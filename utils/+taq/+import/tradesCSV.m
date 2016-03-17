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
    opt.Nrec = 1e5;
    opt.Nblk = 50;
    opt.Scan = {'Delimiter',',','HeaderLines',1};
    opt.Fmt  = 'T%05d';
end
origNblk = opt.Nblk;
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
data   = cell(opt.Nblk,11);
symbol = cell(opt.Nblk,1);
index  = cell(opt.Nblk,1);
ii     = 0;

% LOOP through each file
for f = 1:numel(filenames)
    filename   = unzip(fullfile(path2zip,filenames{f}),tempdir());
    fid        = fopen(filename{end});
    cleanup    = onCleanup(@() finallyCleanup(fid, char(filename)));
    fileIsOpen = true;
    fprintf('%-40s%s\n',filenames{f},datestr(now,'dd HH:MM:SS'))

    while fileIsOpen
        while ii < opt.Nblk && ~feof(fid)
            ii = ii + 1;
            try
                %    1       2       3-5       6       7        8           9           10          11
                % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
                data(ii,:) = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',...
                    opt.Nrec,opt.Scan{:});
                offset     = ftell(fid);
                % Make sure to keep whole day on same mat file
                if ii == opt.Nblk
                    if feof(fid)
                        ressymbol = cell(0,1);
                        residx    = array2table(zeros(0,4), 'VariableNames',{'Id','Date','From','To'});
                        resdata   = cell(1,11);
                    else
                        % Start of last date
                        from = find(diff(data{ii,2}),1,'last')+1;

                        % Day might not have ended yet (whole bulk import has
                        % one date only). Use start of last symbol
                        if isempty(from)
                            [~,~,subs] = unique(data{ii,1});
                            from       = find(diff(subs),1,'last')+1;
                        end

                        % If still empty, import more data
                        if isempty(from)
                            opt.Nblk = opt.Nblk+1;

                            % Defer chunk with last change of date into next bulk import
                            % Note: even if last row is a new date, defer it. We don't
                            %       know if that is the last row for that date.
                        else
                            resdata                    = arrayfun(@(x) data{ii,x}(from:end,:), 1:11,'un',0);
                            data(ii,:)                 = arrayfun(@(x) data{ii,x}(1:from-1,:), 1:11,'un',0);
                            [ressymbol,residx,resdata] = processDataset(resdata);
                        end
                    end
                end
                [symbol{ii},index{ii},data(ii,:)] = processDataset(data(ii,:));

            catch err
                fname = fullfile(path2zip, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
                save(fname,'err')
                rethrow(err)
            end
        end

        % Reset to original block size
        opt.Nblk = origNblk;

        % Close .csv and clean if finished parsing it
        if feof(fid)
            delete(cleanup)
            fileIsOpen = false;
            continue
        end

        % Consolidate parsed blocks into single variables
        [symbol,index,data] = consolidateDataset(symbol,index,data);

        matnum    = matnum+1;
        datafname = fullfile(outdir, sprintf([opt.Fmt '.mat'],matnum));
        idxfname  = fullfile(outdir, sprintf([opt.Fmt '.idx'],matnum));
        save(datafname,'data','-v7.3')
        save(idxfname ,'index','symbol','-v6')
        fprintf('%-40s%s\n',sprintf([opt.Fmt '.mat'],matnum),datestr(now,'dd HH:MM:SS'))

        % Reset containers adding deferred chunks
        data   = [resdata;  cell(opt.Nblk,11)];
        symbol = [{ressymbol}; cell(opt.Nblk,1)];
        index  = [{residx}; cell(opt.Nblk,1)];
        ii     = 1;
    end
end

% LAST iteration exits without saving
[symbol,index,data] = consolidateDataset(symbol,index,data);
matnum              = matnum+1;
datafname           = fullfile(outdir, sprintf([opt.Fmt '.mat'],matnum));
idxfname            = fullfile(outdir, sprintf([opt.Fmt '.idx'],matnum));
save(datafname,'data','-v7.3')
save(idxfname ,'index','symbol','-v6')
fprintf('%-40s%s\n',sprintf([opt.Fmt '.mat'],matnum),datestr(now,'dd HH:MM:SS'))
end

% Process parsed blocks from .csv
function [tck,index,data] = processDataset(data)
nrecords        = size(data{1},1);
% Unique pairs id-date
[tck,~,id]      = unique(data{1});
[iddate,~,subs] = unique([id data{2}],'rows');

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

function finallyCleanup(fid, csvfilename)
fclose(fid);
pause(0.5);
delete(csvfilename);
end