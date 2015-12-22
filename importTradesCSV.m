function matnum = importTradesCSV(path2zip, outdir, isManual, matnum, opt)
% IMPORTTRADESCSV Import zipped .csv of TAQ trades downloaded from WRDS
%
%   IMPORTTRADESCSV(PATH2ZIP, OUTDIR)
%       Imports zipped .CSVs from the PATH2ZIP folder and saves them as
%       .mat and .mst files in OUTDIR.
%
% The .CSVs should be downloaded from WRDS with the following columns. The
% format is on second line (if relevant):
%
%   SYMBOL |   DATE   |   TIME   | PRICE | SIZE | G127 | CORR | COND | EX
%            yyyymmdd   HH:MM:SS
if nargin < 5
    opt.Nrec = 1e5;
    opt.Nblk = 50;
    opt.Scan = {'Delimiter',',','HeaderLines',1};
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
data = cell(opt.Nblk,11);
ids  = cell(opt.Nblk,1);
mst  = cell(opt.Nblk,1);
ii   = 0;

% LOOP through each file
for f = 1:numel(filenames)
    filename   = unzip(fullfile(path2zip,filenames{f}),path2zip);
    fid        = fopen(filename{end});
    cleanup    = onCleanup(@() finallyCleanup(fid, char(filename)));
    fileIsOpen = true;
    fprintf('%-40s%s\n',filenames{f},datestr(now,'dd HH:MM:SS'))

    while fileIsOpen
        while ii < opt.Nblk && ~feof(fid)
            ii     = ii + 1;
            offset = ftell(fid);
            try
                %    1       2       3-5       6       7        8           9           10          11
                % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
                data(ii,:) = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',...
                    opt.Nrec,opt.Scan{:});

                % Make sure to keep whole day on same mat file
                if ii == opt.Nblk
                    if feof(fid)
                        resids = cell(0,1);
                        resmst = array2table(zeros(0,4), 'VariableNames',{'Id','Date','From','To'});
                    else
                        % Start of last date
                        from = find(diff(data{ii,2}),1,'last')+1;

                        % Day might not have ended yet (whole bulk import has
                        % one date only). Use start of last symbol 
                        if isempty(from)
                            [~,~,subs] = unique(data{ii,1});
                            from = find(diff(subs),1,'last')+1;
                        end
                        
                        % If still empty, import more data
                        if isempty(from)
                            opt.Nblk = opt.Nblk+1;

                            % Defer chunk with last change of date into next bulk import
                        elseif from ~= size(data{ii,1},1)
                            resdata                 = arrayfun(@(x) data{ii,x}(from:end,:), 1:11,'un',0);
                            data(ii,:)              = arrayfun(@(x) data{ii,x}(1:from-1,:), 1:11,'un',0);
                            [resids,resmst,resdata] = processDataset(resdata);
                        end
                    end
                end
                [ids{ii},mst{ii},data(ii,:)] = processDataset(data(ii,:));

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
        [ids,mst,data] = consolidateDataset(ids,mst,data);

        matnum    = matnum+1;
        datafname = fullfile(outdir, sprintf('T%04d.mat',matnum));
        mstfname  = fullfile(outdir, sprintf('T%04d.mst',matnum));
        save(datafname,'data','-v7.3')
        save(mstfname ,'mst','ids','-v6')
        fprintf('%-40s%s\n',sprintf('T%04d.mat',matnum),datestr(now,'dd HH:MM:SS'))

        % Reset containers adding deferred chunks
        data = [resdata;  cell(opt.Nblk,11)];
        ids  = [{resids}; cell(opt.Nblk,1)];
        mst  = [{resmst}; cell(opt.Nblk,1)];
        ii   = 1;
    end
end

% LAST iteration exits without saving
[ids,mst,data] = consolidateDataset(ids,mst,data);
matnum         = matnum+1;
datafname      = fullfile(outdir, sprintf('T%04d.mat',matnum));
mstfname       = fullfile(outdir, sprintf('T%04d.mst',matnum));
save(datafname,'data','-v7.3')
save(mstfname ,'mst','ids','-v6')
fprintf('%-40s%s\n',sprintf('T%04d.mat',matnum),datestr(now,'dd HH:MM:SS'))
end

% Process parsed blocks from .csv
function [tck,mst,data] = processDataset(data)
nrecords        = size(data{1},1);
% Unique pairs id-date
[tck,~,id]      = unique(data{1});
[iddate,~,subs] = unique([id data{2}],'rows');

% Starting positions in data
frompos = [0; find(diff(subs))]+1;
topos   = [frompos(2:end)-1; nrecords];

% Store in table and free memory
mst       = table(iddate(:,1),iddate(:,2), frompos, topos, 'VariableNames',{'Id','Date','From','To'});
data(1:2) = {[]};

% Convert to char (and eventually pad with a column of blank)
data{10} = char(data{10});
szCond   = size(data{10});
data{10} = [data{10} repmat(' ',szCond(1),2-szCond(2))];

% Expand tck per date
tck = tck(iddate(:,1));
end

% Consolidate blocks from .csv
function [tck,mst,data] = consolidateDataset(tck,mst,data)
% Count number of records per block
blkcountMst  = cellfun(@(x) size(x,1),mst);
blkcountData = cellfun(@(x) size(x,1),data(:,3));

% Concatenate blocks
data = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
tck  = cat(1,tck{:});
mst  = cat(1,mst{:});

% Some data manipulations
data{1} = cat(2,data{1:3});
data{6} = cat(2,data{6:7});

% Store as table
data = table(data{[1,4:6,8:9]},...
    'VariableNames',{'Time','Price','Volume','G127_Correction','Condition','Exchange'});

% Translate block indexing of From/To to whole file indexing
shift    = RunLength(cumsum([0; blkcountData(1:end-1)]), blkcountMst);
mst.From = mst.From + shift;
mst.To   = mst.To   + shift;

% Re-map Id
[tck,~,mst.Id] = unique(tck);

% Optimize data storage
mst.Id   = uint16(mst.Id);
mst.From = uint32(mst.From);
mst.To   = uint32(mst.To);

% Sort by from
mst = sortrows(mst,'From');

% Glue split series
pos          = find(mst.Id(2:end) == mst.Id(1:end-1) & mst.Date(2:end) == mst.Date(1:end-1));
mst.To(pos)  = mst.To(pos+1);
mst(pos+1,:) = [];
end

function finallyCleanup(fid, csvfilename)
fclose(fid);
pause(0.5);
delete(csvfilename);
end