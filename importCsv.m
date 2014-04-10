%% FINAL files format

% 'data', contains trade data in a 1x6 cell:
%
% 1     [hh, mm, ss]                            uint8
% 2     price                                   single
% 3     size                                    uint32
% 4     [G127 rule,  correction]                uint16
% 5     condition                               2char
% 6     exchange                                1char

% 'mst', master file with ticker symbol, dates, and starting positions:
%
% 1     tickers                                                   cellstr
% 2     [idxTickers, yyyymmdd, startPosData, endPosData, #file]	  uint32


%% Import CSV
dd   = '..\';
%d    = struct('name','7083687f89f5d708.csv');
d    = dir(fullfile(dd,'raw','*.zip'));
opt  = {'Delimiter',',','HeaderLines',1};
N    = 1e5;
name = '';
nite = 200; % approx 400MB in RAM
data = cell(nite,11);
mst  = cell(nite,2);
ii   = 0;
c    = 0;
% d    = d(217:217+35);

% Loop each csv
for f = 1:numel(d)
    filename = unzip(fullfile(dd,'raw',d(f).name),dd);
    %filename = {'E:\data\7083687f89f5d708.csv'};
    fid      = fopen(filename{end});
    status   = true;
    fprintf('%-40s%s\n',d(f).name,datestr(now,'dd HH:MM:SS'))
    
    while status
        while ii < nite && ~feof(fid)
            ii = ii + 1;
            %    1       2       3-5       6       7        8           9           10          11
            % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
            data(ii,:)         = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',N,opt{:});
            % Unique tickers
            [mst{ii,1},~,idx]  = unique(data{ii,1});
            % Unique pairs ticker-date
            [symbdate,~,stPos] = unique([idx data{ii,2}],'rows');
            mst{ii,2}          = [symbdate, find(diff([0;stPos]))];
            % Clear first two columns for memory
            data(ii,1:2)       = {[]};
            % Convert to char
            data{ii,10}        = char(data{ii,10});
            szCond             = size(data{ii,10});
            data{ii,10}        = [data{ii,10} repmat(' ',szCond(1),2-szCond(2))];
            % Add end positions to the pair id-data before start positions get shuffled in the re-indexing
            mst{ii,2}(:,end+1) = mst{ii,2}(:,end) + diff([mst{ii,2}(:,end);szCond(1)+1]) - 1; 
        end
        
        % Close and clean
        if feof(fid) && ii < 200
            status = fclose(fid);
            pause(0.5)
            delete(filename{:})
            continue
        end
        
        % Collect everything in one cell per column data and reorganize a bit
        lenData           = uint32(cellfun('size',data(:,3),1));
        lenMst            = uint32(cellfun('size',mst,1));
        data              = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
        data{1}           = cat(2,data{1:3});
        data{6}           = cat(2,data{6:7});
        data([2:3 7])     = [];
        mst               = {cat(1,mst{:,1}) cat(1,mst{:,2})};
        
        % Translate block indexing to whole file indexing
        idx               = cumsum(lenMst(1:end-1,2))+1;                % Start position of each block
        val               = [lenData(1:ii-1) - mst{2}(idx-1,3)+1,... 
                             mst{2}(idx,end)];                          % Calculate within blocks last-position increments
        mst{2}(:,[1,3:4]) = [1,1,mst{2}(1,4); diff(mst{2}(:,[1,3:4]))]; % Diff positions to get increments
        mst{2}(idx,1)     = 1;                                          % Reset 0s from diff (recall it's uint type) to 1, i.e. start pos to 1
        mst{2}(idx,3:4)   = val;                                        % Set all last increments within block
        mst{2}(:,[1 3:4]) = cumsum(mst{2}(:,[1,3:4]));                  % Cumulate back to get file positions
        
        % Remove ticker repetitions and redundant ticker-date pairs retaining only first start pos and max end pos
        % pos:             tickerNew -> tickerOld (idx)        i.e. maps unique tickers to positions with repetitions s.t. tickerNew(pos) = tickerOld
        % mst{2}(:,1):     tickerOld -> id (idx)               i.e. maps repeated tickers to positions (simple enumeration)
        % pos(mst{2}(:,1)) tickerOld -> id (idx) -> tickerNew  i.e. different positions for repeated tickers are mapped to unique position per ticker
        [mst{1},~,posGrow]    = unique(mst{1});
        mst{2}(:,1)       = posGrow(mst{2}(:,1));
        [~,idx,subs]      = unique(mst{2}(:,1:2),'rows','first');
        mst{2}            = [mst{2}(idx,1:3) accumarray(subs,mst{2}(:,4),[],@max)];
        
        % Saving cell data and master as is - tested 1.6GB - 21 sec to save; 4.5 sec to load; -30% space
        c = c+1;
        save(fullfile(dd,sprintf('T%04d.mat',c)),'data','mst','-v7.3')
        fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
        % Reset containers
        data = cell(nite,11);
        mst  = cell(nite,2);
        ii   = 0;
    end
end

% LAST iteration exits without processing saving, thus do it here
lenData           = uint32(cellfun('size',data(1:ii,3),1));
lenMst            = uint32(cellfun('size',mst(1:ii,:),1));
data              = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
data{1}           = cat(2,data{1:3});
data{6}           = cat(2,data{6:7});
data([2:3 7])     = [];
mst               = {cat(1,mst{:,1}) cat(1,mst{:,2})};
idx               = cumsum(lenMst(1:end-1,2))+1;
val               = [lenData(1:end-1)-mst{2}(idx-1,3)+1, mst{2}(idx,end)];
mst{2}(:,[1,3:4]) = [1, 1, mst{2}(1,4); diff(mst{2}(:,[1, 3:4]))];
mst{2}(idx,1)     = 1;
mst{2}(idx,3:4)   = val;
mst{2}(:,[1 3:4]) = cumsum(mst{2}(:,[1,3:4]));
[mst{1},~,posGrow]    = unique(mst{1});
mst{2}(:,1)       = posGrow(mst{2}(:,1));
[~,idx,subs]      = unique(mst{2}(:,1:2),'rows','first');
mst{2}            = [mst{2}(idx,1:3) accumarray(subs,mst{2}(:,4),[],@max)];
c = c+1;
save(fullfile(dd,sprintf('T%04d.mat',c)),'data','mst','-v7.3')
fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
% Set data to read only and hidden
fileattrib(fullfile(dd,'*.mat'),'+h -w','','s')

%% Generate grand master.mat file
% Read .mat filenames
dd = '..\';
d  = dir(fullfile(dd,'*.mat'));

tic
% Open matlabpool
if matlabpool('size') == 0
    matlabpool 4
end
toc
tic
% Preallocate and load mst
nfiles = numel(d);
mst    = cell(nfiles,2);
for f = 1:nfiles
    s = load(fullfile(dd,d(f).name));
    mst(f,:)  = s.mst;
end
matlabpool close
toc
% Retrieve number of rows per block
len            = cellfun('size',mst(:,2),1);
% Add number of file as 4th column to mst{2} and cat all mst
mst(:,2)       = arrayfun(@(x) [mst{x,2} repmat(x,len(x),1)],1:numel(len),'un',0);
mst            = {cat(1,mst{:,1}), cat(1,mst{:,2})}; 
% Re-index id to whole block
posGrow            = cumsum(len(1:end-1));
mst{2}(:,1)    = [1; diff(mst{2}(:,1))];
mst{2}(posGrow+1)  = 1;
mst{2}(:,1)    = cumsum(mst{2}(:,1));
% Make unique tickers
[mst{1},~,posGrow] = unique(mst{1});
mst{2}(:,1)    = posGrow(mst{2}(:,1));
% Sort according to id-date pair
mst{2}         = sortrows(mst{2},1:2);
% Save
save(fullfile(dd,'master.mat'),'mst','-v6')
%% Check if imported correctly
tic;load('..\master','-mat');toc
for y = 1993:2010
    dd = sprintf('..\\datachk\\%d\\',y);
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
        
        idxM = mst{2}(:,1) == find(strcmpi(regexp(d(f).name,'\w+(?=_)','match'), mst{1})) & ...
               fix(mst{2}(:,2)/1e4) == y;
        tmp  = mst{2}(idxM,:);
        cds  = unique(tmp(:,end));
        nCds = numel(cds);
        dataC = cell(nCds,6);
        for ii = 1:nCds
            s = load(fullfile('..\',sprintf('T%04d.mat',cds(ii))));
            dataC(ii,:) = s.data;
            idxTmp = tmp(:,end) == cds(ii);
            idx    = false(size(dataC{ii,2}));
            for it = find(idxTmp)'
                idx(tmp(it,3):tmp(it,4)) = true;
            end
            dataC(ii,:) = arrayfun(@(x) dataC{ii,x}(idx,:),1:6,'un',0);
        end
        dataC = arrayfun(@(x) cat(1,dataC{:,x}),1:6,'un',0);
        disp(isequal(dataC,data(3:end)))
    end
end
%% Convert into dataset (added later)
% Convert overall master file
d   = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')
ids = mst{1};
mst = dataset({uint16(mst{2}(:,1)),'Id'},...
              {mst{2}(:,2:4),'Date','From','To'},...
              {uint16(mst{2}(:,5)) ,'File'});
save(fullfile(d,'master.mat'),'mst','ids','-v6')
movefile(fullfile(d,'master.mat'),fullfile(d,'master'))

% Convert all data and intra file master
dd = dir(fullfile(d,'*.mat'));
for f = 1:numel(dd)
    tic
    fname = fullfile(d,dd(f).name);
    disp(fname)
    s     = load(fname);
    ids   = s.mst{1};
    mst   = dataset({uint16(s.mst{2}(:,1)),'Id'}, {s.mst{2}(:,2:4),'Date','From','To'});
    % Cell2dataset is slower
    data  = dataset({s.data{1},'Time'},...
                    {s.data{2},'Price'},...
                    {s.data{3},'Volume'},...
                    {s.data{4},'G127_Correction'},...
                    {s.data{5},'Condition'},...
                    {s.data{6},'Exchange'});
    % This ensures that the original is not overwritten partially in case of error
    datname = regexprep(fname,'\.mat','\.dat');
    save(datname,'data','mst','ids')
    toc
    delete(fname)
    movefile(datname,fname)
    toc
end
fileattrib(fullfile(d,'*.mat'),'+h -w','','s')
%% Avoid splits
tic
% List of files
cd C:\HFbetas
d  = '.\data\TAQ';
dd = dir(fullfile(d,'*.mat'));

% Enable write
fileattrib(fullfile(d,'*.mat'),'-h +w','','s')

% Find splits: a daily series can be split among 2 files (not more, checked)
load(fullfile(d,'master'),'-mat')
[unSer,~,subs] = unique(mst(:,{'Id','Date'}));

% Check if the split spans more than 2 files [OK - it does not]
if nnz(accumarray(subs,1) > 2); error('Split spans more than 2 files'); end

splitpos   = sort(find(diff(subs) == 0));
grow       = [dataset({ids(mst.Id(splitpos)),'Name'}), mst(splitpos,:)];
shorten    = [dataset({ids(mst.Id(splitpos+1)),'Name'}), mst(splitpos+1,:)];

% Check the split occurs at end and beginnig of the 2 contiguous files [OK - always the case]
maxfiles = accumarray(mst.File, mst.To,[],@max);
if nnz(grow.To ~= maxfiles(mst.File(splitpos))) || nnz(shorten.From ~= 1) 
    error('Split does not happen at the end/beginning of two contiguous files')
end

% LOOP (no-brainer)
for ii = size(grow,1):-1:1
    % Load the two files among which a daily series is split
    filenameA = fullfile(d, sprintf('T%04d.mat', grow.File(ii)));
    A = load(filenameA);
    filenameB = fullfile(d, sprintf('T%04d.mat', shorten.File(ii)));
    B = load(filenameB);
    
    % Append to A (note idx always starts at 1, checked above)
    idx    = 1:shorten.To(ii);
    nidx   = numel(idx);
    A.data = [A.data; B.data(idx,:)];
    
    % Edit data B
    B.data(idx,:) = [];
    
    % Edit mst of A (grow only the To attribute of the split series)
    posGrow      = find(ismember([dataset({A.ids(A.mst.Id),'Name'}), A.mst(:,'Date')],...
                                  grow(ii,{'Name','Date'})));
    [~,checkPos] = max(A.mst.To);
    if posGrow ~= checkPos; error('Wrong position to grow A');  end
    A.mst.To(posGrow) = A.mst.To(posGrow) + nidx;
    
    % Edit mst of B (shift everything by the length of the split series)
    posShorten      = find(ismember([dataset({B.ids(B.mst.Id),'Name'}), B.mst(:,'Date')],...
                                     shorten(ii,{'Name','Date'})));
    [~,checkPos] = min(B.mst.From);
    if posShorten ~= checkPos; error('Wrong position to shorten B');  end
    shortId             = B.mst.Id(posShorten);
    B.mst(posShorten,:) = [];
    B.mst.From          = B.mst.From - nidx;
    B.mst.To            = B.mst.To   - nidx;
    
    % Edit ids of B (if the split series does not have any other days, disappears)
    if ~any(B.mst.Id == shortId)
        B.mst.Id(B.mst.Id > shortId) = B.mst.Id(B.mst.Id > shortId) -1;
        B.ids(shortId) = [];
    end
    
    % Edit master
    mst.To(splitpos(ii),:) = A.mst.To(posGrow);
    mst(splitpos(ii)+1,:)  = [];
    ifile                  = mst.File == shorten.File(ii);
    mst.From(ifile)        = mst.From(ifile) - nidx;
    mst.To(ifile)          = mst.To(ifile) - nidx;
    
    % Check that changes are coherent
    if ~isequal(ids(mst.Id(ifile)), B.ids(B.mst.Id));   error('B: problem with IDs.');     end
    if ~isequal(mst(ifile,2:end-1), B.mst(:,2:end));    error('B: problem with mst.');    end
    ifile = mst.File == grow.File(ii);
    if ~isequal(ids(mst.Id(ifile)), A.ids(A.mst.Id));   error('A: prblem with IDs.');     end
    if ~isequal(mst(ifile,2:end-1), A.mst(:,2:end));    error('A: problem with mst.');    end
    
    % Save progress
    save(filenameA,'-struct','A','ids','mst','data')
    save(filenameB,'-struct','B','ids','mst','data')
    toc
end
save(fullfile(d,'master'),'mst','ids','-mat','-v6')
fileattrib(fullfile(d,'*.mat'),'+h -w','','s')
sendolmail('o.komarov11@imperial.ac.uk',sprintf('Terminated in %.0f',toc),'')
%% Check grand master file against checkdata
cd C:\HFbetas
localpath = '.\data\TAQ';
load(fullfile(localpath,'master'),'-mat')
externalpath = 'D:\HFbetas\data\TAQ\datachk';
for y = 1993:2010
    dd = fullfile(externalpath,sprintf('%d',y));
    d  = dir(fullfile(dd,'*.zip'));
    disp(dd)
    for f = 1:numel(d)
        filename   = unzip(fullfile(dd,d(f).name),dd);
        fid        = fopen(filename{end});
        data       = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c','Delimiter',',','HeaderLines',1);
        data{1,10} = char(data{1,10});
        szCond     = size(data{1,10});
        data{1,10} = [data{1,10} repmat(' ',szCond(1),2-szCond(2))];
        data       = dataset({cat(2,data{3:5}), 'Time'}, {data{6},'Price'},{data{7},'Volume'},...
                             {cat(2,data{8:9}),'G127_Correction'},{data{10},'Condition'},{data{11},'Exchange'});
        fclose(fid);
        delete(filename{1});
        
        idxM = mst.Id == find(strcmpi(regexp(d(f).name,'\w+(?=_)','match'), ids)) & ...
               fix(mst.Date/1e4) == y;
        tmp  = mst(idxM,:);
        cds  = unique(tmp.File);
        nCds = numel(cds);
        dataC = cell(nCds,1);
        for ii = 1:nCds
            s = load(fullfile(localpath,sprintf('T%04d.mat',cds(ii))));
            dataC{ii} = s.data;
            idxTmp = tmp.File == cds(ii);
            idx    = false(size(dataC{ii}.Price));
            for it = find(idxTmp)'
                idx(tmp.From(it):tmp.To(it)) = true;
            end
            dataC{ii} = dataC{ii}(idx,:);
        end
        dataC = cat(1,dataC{:});
        disp(isequal(dataC,data))
    end
end

% Now check that single master files compose into grand
% Read .mat filenames
d  = dir(fullfile(localpath,'*.mat'));

% Open matlabpool
if matlabpool('size') == 0
    matlabpool 4
end

% Preallocate and load mst
nfiles    = numel(d);
[mst2,ids2] = deal(cell(nfiles,1));
tic
for f = 1:nfiles
    disp(f)
    s = load(fullfile(localpath,d(f).name));
    mst2{f} = s.mst;
    ids2{f} = s.ids;
end
matlabpool close
toc
% Retrieve number of rows per block
len               = cellfun(@(x)size(x,1),mst2);
% Add number of file as 4th column to mst{2} and cat all mst
ids2               = cat(1,ids2{:});
mst2               = cat(1,mst2{:});
mst2.Id            = uint32(mst2.Id);
mst2.File          = RunLength((1:numel(len))', len);
% Re-index id to whole block
posGrow           = cumsum(len(1:end-1));
mst2.Id            = [1; diff(mst2.Id)];
mst2.Id(posGrow+1) = 1;
mst2.Id            = cumsum(mst2.Id);
% Make unique tickers
[ids2,~,posGrow]   = unique(ids2);
mst2.Id            = posGrow(mst2.Id);
% Sort according to id-date pair
mst2               = sortrows(mst2,{'Id','Date'});
% Compare
load(fullfile(localpath,'master'),'-mat')
isequal(mst2,mst), isequal(ids2,ids)
%% Sort mst with From
cd C:\HFbetas
localpath = '.\data\TAQ';

d  = dir(fullfile(localpath,'*.mat'));

% Open matlabpool
if matlabpool('size') == 0
    matlabpool 4
end

% Preallocate, load mst and sort
nfiles    = numel(d);
tmp = cell(nfiles,1);
tic
parfor f = 1:nfiles
    disp(f)
    filename = fullfile(localpath,d(f).name);
    s = load(filename,'mst');
    tmp{f} = sortrows(s.mst,'From');
end
toc
% Overwrite mst with append or I lose the other vars
fileattrib(fullfile(localpath,'*.mat'),'-h +w','','s')
for f = 1:nfiles
    disp(f)
    mst = tmp{f};
    save(fullfile(localpath,d(f).name),'mst','-append')
end
fileattrib(fullfile(localpath,'*.mat'),'+h -w','','s')