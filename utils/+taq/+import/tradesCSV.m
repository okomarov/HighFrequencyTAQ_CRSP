function matnum = tradesCSV(path2zip, outdir, isManual, matnum, opt)
% TRADESCSV Import zipped .csv of TAQ trades downloaded from WRDS
%
%   TRADESCSV(PATH2ZIP, OUTDIR)
%       Imports zipped .CSVs from the PATH2ZIP folder and saves them as
%       .mat and .idx files in OUTDIR.
%
%   Optional inputs:
%       ISMANUAL - boolean (default: false). Set to true if the .CSVs were
%                  downloaded from WRDS in the format specified below:
%
%           symbol |   date   |   time   | price | size | g127 | corr | cond | ex
%                    yyyymmdd   HH:MM:SS
%
%       MATNUM - numeric scalar (default: 0). Starts the numbering of
%                .mat/idx files from MATNUM + 1.
%
%       OPT - structure with options (defaults: see in the code).
%           .Nrec - number of records with each textscan read.
%           .Nblk - number of reads before saving results to .mat/idx
%           .Scan - textscan options
%           .Fmt  - read format for textscan
%
% See also: TEXTSCAN,TRADESDVD

if nargin < 5
    opt.Nrec        = 1e5;
    opt.Nblk        = 50;
    opt.Scan        = {'Delimiter',',','HeaderLines',1};
    opt.FilenameFmt = 'T%05d';
end
if nargin < 4 || isempty(matnum)
    matnum = 0;
end

% Read list of .zip files
d = dir(fullfile(path2zip,'*.zip'));
if isManual
    filenames = regexp(sort({d.name})', '\d{6}_\d{1,2}.zip','match','once');
else
    filenames = regexp(sort({d.name})', 'CT_\d{8}.zip','match','once');
end
filenames = filenames(~cellfun('isempty', filenames));

% Preallocate
[data,symbol,index, resdata] = preallocateContainters(opt);
ii                           = 0;

% LOOP through each file
for f = 1:numel(filenames)
    filename = unzip(fullfile(path2zip,filenames{f}),tempdir());
    fid      = fopen(filename{end});
    cleanup  = onCleanup(@() finallyCleanup(fid, char(filename)));
    fprintf('%-40s%s\n',filenames{f},datestr(now,'dd HH:MM:SS'))

    while ~feof(fid)
        while ii < opt.Nblk
            try
                ii = ii + 1;
                [data(ii,:), symbol{ii}, index{ii}, resdata, opt] = importBlock(fid,ii,opt,resdata);

                unzipNewFile = feof(fid) && ii < fix(opt.Nblk * 0.8);
                if unzipNewFile
                    break
                elseif feof(fid)
                    % File is over and number of bulk import not complete
                    % BUT we consolidate anywas to set a breakpoint in the
                    % importing routine
                    fprintf('%-40s BREAKPOINT\n',filenames{f})
                    ii = opt.Nblk;
                end
            catch err
                fname = fullfile(path2zip, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
                save(fname,'err')
                rethrow(err)
            end
        end

        if unzipNewFile
            break
        end

        [symbol,index,data] = consolidateDataset(symbol,index,data);

        matnum = saveFiles(outdir,matnum,opt,symbol,index,data);

        [data,symbol,index] = preallocateContainters(opt);
        ii                  = 0;
    end
end
if ~isempty(data)
    [symbol,index,data] = consolidateDataset(symbol,index,data);
    matnum = saveFiles(outdir,matnum,opt,symbol,index,data);
end
end

% Imports block of data from csv
function [data, symbol, index, resdata, opt] = importBlock(fid,iter,opt,resdata)

%    1       2       3-5       6       7        8           9           10          11
% ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
fmt  = '%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c';
data = textscan(fid, fmt, opt.Nrec,opt.Scan{:});

% For debugging
opt.Offset = ftell(fid);

% Preallocate empty residual
if ~isempty(resdata)
    data            = cellfun(@(res,d) [res;d], resdata,data,'un',0);
    [~,~,~,resdata] = preallocateContainters(opt);
end

% Determine residual data that goes into next bulk import
if iter >= opt.Nblk

    % Ensure the whole day is within one mat file
    % Start of last date
    from = find(diff(data{2}),1,'last')+1;

    % Day might not have ended yet (whole bulk import has
    % one date only). Use start of last symbol
    if isempty(from)
        [~,~,subs] = unique(data{1});
        from       = find(diff(subs),1,'last')+1;
    end

    % If still empty, import more data
    if isempty(from)
        [data, symbol, index, resdata, opt] = importBlock(fid,iter+1,opt, data);
        return
        % Defer chunk with last change of date into next bulk import
        % Note: even if last row is a new date, defer it. We don't
        %       know if that is the last row for that date.
    else
        resdata = arrayfun(@(x) data{x}(from:end,:), 1:11,'un',0);
        data    = arrayfun(@(x) data{x}(1:from-1,:), 1:11,'un',0);
    end
end
[symbol,index,data] = processDataset(data);
end

% Process parsed blocks from .csv
function [tck,index,data] = processDataset(data)
nrecords        = size(data{1},1);
% Unique pairs id-date
[tck,~,id]      = unique(data{1});
[iddate,~,subs] = unique([id data{2}],'rows','stable');

% Starting positions in data
frompos = [0; find(diff(subs))]+1;
topos   = [frompos(2:end)-1; nrecords];

% Store in table and free memory
index     = table(iddate(:,1),iddate(:,2), frompos, topos, 'VariableNames',{'Id','Date','From','To'});
data(1:2) = {[]};

% Convert to char (and eventually pad with a column of blank)
data{10} = char(data{10});
szCond   = size(data{10});
data{10} = [data{10} repmat(' ',szCond(1),2-szCond(2))];

% Expand tck per date
tck = tck(iddate(:,1));
end

% Consolidate blocks from .csv
function [tck,index,data] = consolidateDataset(tck,index,data)
% Count number of records per block
blkcountMst  = cellfun(@(x) size(x,1),index);
blkcountData = cellfun(@(x) size(x,1),data(:,3));

% Concatenate blocks
data  = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
tck   = cat(1,tck{:});
index = cat(1,index{:});

% Some data manipulations
data{1} = cat(2,data{1:3});
data{6} = cat(2,data{6:7});

% Store as table
data = table(data{[1,4:6,8:9]},...
    'VariableNames',{'Time','Price','Volume','G127_Correction','Condition','Exchange'});

% Translate block indexing of From/To to whole file indexing
shift      = RunLength(cumsum([0; blkcountData(1:end-1)]), blkcountMst);
index.From = index.From + shift;
index.To   = index.To   + shift;

% Re-map Id
[tck,~,index.Id] = unique(tck);

% Optimize data storage
index.Id   = uint16(index.Id);
index.From = uint32(index.From);
index.To   = uint32(index.To);

% Sort by from
index = sortrows(index,'From');

% Glue split series
key           = uint64(index.Id)*1e8 + uint64(index.Date);
[~,pos,subs]  = unique(key,'first');
to            = accumarray(subs, index.To,[],@max);
index.To(pos) = to;
index         = index(pos,:);
index         = sortrows(index,'From');
end

% Saves .mat and .idx
function matnum = saveFiles(outdir, matnum, opt, symbol, index, data)
matnum = matnum+1;

datafname = fullfile(outdir, sprintf([opt.FilenameFmt '.mat'],matnum));
save(datafname,'data','-v7.3')

idxfname = fullfile(outdir, sprintf([opt.FilenameFmt '.idx'],matnum));
save(idxfname ,'index','symbol','-v6')

fprintf('%-40s%s\n',sprintf([opt.FilenameFmt '.mat'],matnum),datestr(now,'dd HH:MM:SS'))
end

% Preallocations
function [data,symb,idx,rdata] = preallocateContainters(opt)
data  = cell(opt.Nblk,11);
symb  = cell(opt.Nblk,1);
idx   = cell(opt.Nblk,1);
rdata = cell(0,11);
end

% Fid and unzipped cleanup
function finallyCleanup(fid, csvfilename)
fclose(fid);
pause(0.5);
delete(csvfilename);
end