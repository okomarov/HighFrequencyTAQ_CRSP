function importCsv2(path2zip,opt)
% IMPORTCSV Imports WRDS .csv with TAQ prices
%   
%   IMPORTCSV(PATH2ZIP) Imports zipped .CSVs from the PATH2ZIP directory
%                       and saves them into .mat files under
%                       $PATH2ZIP$\mat\T####.mat. 
%
%
% The .CSVs should be downloaded from WRDS with the following columns. The 
% format is on second line (if relevant):
%
%   SYMBOL |   DATE   |   TIME   | PRICE | SIZE | G127 | CORR | COND | EX
%            yyyymmdd   HH:MM:SS
%                       
% The .mat files will contain the variables data (table), mst (table) and
% ids (cellstring).
%
% - data:
%       .Time [hh, mm, ss]                          uint8
%       .Price                                      single
%       .Volume                                     uint32
%       .G127_Correction [G127 rule,  correction]   uint16
%       .Condition                                  2char
%       .Exchange                                   1char
%
% - ids, unique list of tickers                     cellstring 
%
% - mst, master of id-date pairs mapping to data and ids:
%       .Id, map to ids                             uint16
%       .Date, yyyymmdd                             uint32
%       .From, starting row in data                 uint32
%       .To,   ending row in data                   uint32
%
% For details on the data and its fields see <a href="http://www.nyxdata.com/Data-Products/Monthly-TAQ-DVD">Monthly TAQ DVD</a>.

if nargin < 2
    opt.Nrec = 1e5;
    opt.Nblk = 100;
    opt.Scan = {'Delimiter',',','HeaderLines',1}; 
end
    
%% Import CSV
d         = dir(fullfile(path2zip,'*.zip'));
filenames = regexp(sort({d.name})', '\d{6}_\d{1,2}.zip','match','once');
filenames = filenames(~cellfun('isempty', filenames));

% Preallocate
data = cell(opt.Nblk,11);
ids  = cell(opt.Nblk,1);
mst  = cell(opt.Nblk,1);
ii   = 0;
c    = 0;

% Loop each csv
for f = 1:numel(filenames)
    filename = unzip(fullfile(path2zip,filenames{f}),path2zip);
    fid      = fopen(filename{end});
    status   = true;
    fprintf('%-40s%s\n',filenames{f},datestr(now,'dd HH:MM:SS'))
    
    while status
        while ii < opt.Nblk && ~feof(fid)
            ii = ii + 1;
            
            %    1       2       3-5       6       7        8           9           10          11
            % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
            data(ii,:) = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',...
                                  opt.Nrec,opt.Scan{:});
            
            % Make sure to keep whole day on same file
            if ii == opt.Nblk
                from       = find(diff(data{ii,2}),1,'last')+1;
                resdata    = arrayfun(@(x) data{ii,x}(from:end,:), 1:11,'un',0);
                data(ii,:) = arrayfun(@(x) data{ii,x}(1:from-1,:), 1:11,'un',0);
                [resids,resmst,resdata] = processDataset(resdata);
            end
            [ids{ii},mst{ii},data(ii,:)] = processDataset(data(ii,:));
        end
        
        % Close .csv and clean if finished parsing it
        if feof(fid)
            status = fclose(fid);
            pause(0.5)
            delete(filename{:})
            continue
        end
        
        % Consolidate parsed blocks into single variables
        [ids,mst,data] = consolidateDataset(ids,mst,data);
        
        % Saving cell data and master as is - tested 1.6GB - 21 sec to save; 4.5 sec to load; -30% space
        c = c+1;
        save(fullfile(path2zip,'mat', sprintf('T%04d.mat',c)),'data','mst','ids','-v7.3')
        fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
        % Reset containers
        data = [resdata;  cell(opt.Nblk,11)];
        ids  = [{resids}; cell(opt.Nblk,1)];
        mst  = [{resmst}; cell(opt.Nblk,1)];
        ii   = 1;
    end
end

% LAST iteration exits without processing saving, thus do it here
[ids,mst,data] = consolidateDataset(ids,mst,data);
c = c+1;
save(fullfile(path2zip,'mat', sprintf('T%04d.mat',c)),'data','mst','ids','-v7.3')
fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
%% Set data to read only and hidden
fileattrib(fullfile(path2zip,'mat','*.mat'),'+h -w','','s')

%% Generate grand master.mat file
makeMasterFile(fullfile(path2zip,'mat'))

%% Check if imported correctly
load(fullfile(path2zip,'mat','master'),'-mat');
for y = 1993:2010
    dd = fullfile(path2zip, '..\datachk\',sprintf('%d',y));
    d  = dir(fullfile(dd,'*.zip'));
    disp(dd)
    for f = 1:numel(d)
        filename   = unzip(fullfile(dd,d(f).name),dd);
        fid        = fopen(filename{end});
        data       = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c','Delimiter',',','HeaderLines',1);
        data{1,10} = char(data{1,10});
        szCond     = size(data{1,10});
        data{1,10} = [data{1,10} repmat(' ',szCond(1),2-szCond(2))];
        data       = [data(1:2) cat(2,data{3:5}) data(6:7) cat(2,data{8:9}) data(10:11)];
        fclose(fid);
        delete(filename{1});
        
        idxM = mst.Id == find(strcmpi(regexp(d(f).name,'\w+(?=_)','match'), ids)) & ...
               fix(mst.Date/1e4) == y;
        tmp  = mst(idxM,:);
        cds  = unique(tmp(:,end));
        nCds = numel(cds);
        dataC = cell(nCds,6);
        for ii = 1:nCds
            s           = load(fullfile(path2zip,'mat',sprintf('T%04d.mat',cds(ii))));
            dataC(ii,:) = s.data;
            idxTmp      = tmp(:,end) == cds(ii);
            idx         = false(size(dataC{ii,2}));
            for it = find(idxTmp)'
                idx(tmp(it,3):tmp(it,4)) = true;
            end
            dataC(ii,:) = arrayfun(@(x) dataC{ii,x}(idx,:),1:6,'un',0);
        end
        dataC = arrayfun(@(x) cat(1,dataC{:,x}),1:6,'un',0);
        disp(isequal(dataC,data(3:end)))
    end
end

%% Check for splits
load(fullfile(path2zip,'mat','master'),'-mat')
[unSer,~,subs] = unique(mst(:,{'Id','Date'}));

% Check if the split spans more than 2 files [OK - it does not]
if nnz(accumarray(subs,1) > 1); error('There are splits'); end
end

% Process parsed blocks from .csv
function [tck,mst,data] = processDataset(data)
% Number of records parsed < nrec if file ends
nrecords        = size(data{1},1);
% Unique tickers
[tck,~,id]      = unique(data{1});
% Unique pairs id-date
[iddate,~,subs] = unique([id data{2}],'rows');
% Starting positions in data
frompos         = [0; find(diff(subs))]+1;
mst             = table(iddate(:,1),iddate(:,2), frompos, 'VariableNames',{'Id','Date','From'});
% End positions (before start positions get shuffled)
mst.To          = [mst.From(2:end)-1; nrecords];
% Clear first two columns for memory
data(1:2)       = {[]};
% Convert to char (and eventually pad with a column of blank)
data{10}        = char(data{10});
szCond          = size(data{10});
data{10}        = [data{10} repmat(' ',szCond(1),2-szCond(2))];
% Expand tck per date
tck             = tck(iddate(:,1));
end

% Consolidate blocks from .csv
function [tck,mst,data] = consolidateDataset(tck,mst,data)
% Count number of records per block
blkcountMst   = cellfun(@(x) size(x,1),mst);
blkcountData  = cellfun(@(x) size(x,1),data(:,3));
% Concatenate blocks
data          = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
tck           = cat(1,tck{:});
mst           = cat(1,mst{:});
% Some data manipulations
data{1}       = cat(2,data{1:3});
data{6}       = cat(2,data{6:7});
% Store as table
data          = table(data{[1,4:6,8:9]},...
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