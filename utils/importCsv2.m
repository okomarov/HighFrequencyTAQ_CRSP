function importCsv2(path2zip)
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
%d    = struct('name','7083687f89f5d708.csv');
d    = dir(fullfile(path2zip,'*.zip'));
opt  = {'Delimiter',',','HeaderLines',1};
nrec = 1e5;
N    = 100; % approx 200MB in RAM

data = cell(N,11);
[tck, mst] = deal(cell(N,1));
ii   = 0;
c    = 0;
% d    = d(217:217+35);

% Loop each csv
for f = 1:numel(d)
    filename = unzip(fullfile(path2zip,d(f).name),path2zip);
    %filename = {'E:\data\7083687f89f5d708.csv'};
    fid      = fopen(filename{end});
    status   = true;
    fprintf('%-40s%s\n',d(f).name,datestr(now,'dd HH:MM:SS'))
    
    while status
        while ii < N && ~feof(fid)
            ii = ii + 1;
            
            %    1       2       3-5       6       7        8           9           10          11
            % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
            data(ii,:) = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',nrec,opt{:});
            
           % Make sure to keep whole day on same file
            if ii == N-1
                from    = find(diff(data{ii,2}),1,'last')+1;
                resdata = arrayfun(@(x) data{ii,x}(from:end,:), 1:11,'un',0);
                data    = arrayfun(@(x) data{ii,x}(1:from-1,:), 1:11,'un',0);
            end
            
            nrecords           = size(data{ii,1},1);
            % Unique tickers
            [tck{ii},~,id]     = unique(data{ii,1});
            % Unique pairs id-date
            [iddate,~,subs]    = unique([id data{ii,2}],'rows');
            % Starting positions in data
            frompos            = [1; find(diff(subs))];
            mst{ii}            = table(iddate(:,1),iddate(:,2), frompos, 'VariableNames',{'Id','Date','From'});
            % End positions (before start positions get shuffled)
            mst{ii}.To         = [mst{ii,2}.From(2:end)-1; nrecords]; 
            % Clear first two columns for memory
            data(ii,1:2)       = {[]};
            % Convert to char 
            data{ii,10}        = char(data{ii,10});
            szCond             = size(data{ii,10});
            data{ii,10}        = [data{ii,10} repmat(' ',szCond(1),2-szCond(2))];
        end
        
        % Close and clean
        if feof(fid) && ii < N
            status = fclose(fid);
            pause(0.5)
            delete(filename{:})
            continue
        end
        
        % Count number of records per block 
        blkcountMst = uint32(cellfun('size',mst,1));
        
        % Concatenate blocks
        data          = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
        tck           = cat(1,tck{:});
        mst           = cat(1,mst{:});
        % Some data manipulations
        data{1}       = cat(2,data{1:3});
        data{6}       = cat(2,data{6:7});
        data([2:3 7]) = [];
        
        % Translate block indexing of From/To to whole file indexing
        shift    = RunLength(cumsum([0; blkcountMst]), blkcountMst);
        mst.From = mst.From + shift;
        mst.To   = mst.To   + shift;
        
        % Re-map Id
        [tck,~,mst.Id] = unique(tck);
        
        % Saving cell data and master as is - tested 1.6GB - 21 sec to save; 4.5 sec to load; -30% space
        c = c+1;
        save(fullfile(path2zip,'mat', sprintf('T%04d.mat',c)),'data','mst','tck','-v7.3')
        fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
        % Reset containers
        data       = [resdata; cell(N,11)];
        [tck, mst] = deal(cell(N,1));
        ii         = 1;
    end
end

%%  Rewrite from here
% LAST iteration exits without processing saving, thus do it here
blkcountData           = uint32(cellfun('size',data(1:ii,3),1));
blkcountMst            = uint32(cellfun('size',mst(1:ii,:),1));
data              = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
data{1}           = cat(2,data{1:3});
data{6}           = cat(2,data{6:7});
data([2:3 7])     = [];
mst               = {cat(1,mst{:,1}) cat(1,mst{:,2})};
idx               = cumsum(blkcountMst(1:end-1,2))+1;
val               = [blkcountData(1:end-1)-mst{2}(idx-1,3)+1, mst{2}(idx,end)];
mst{2}(:,[1,3:4]) = [1, 1, mst{2}(1,4); diff(mst{2}(:,[1, 3:4]))];
mst{2}(idx,1)     = 1;
mst{2}(idx,3:4)   = val;
mst{2}(:,[1 3:4]) = cumsum(mst{2}(:,[1,3:4]));
[tck,~,posGrow]    = unique(tck);
mst{2}(:,1)       = posGrow(mst{2}(:,1));
[~,idx,subs]      = unique(mst{2}(:,1:2),'rows','first');
mst{2}            = [mst{2}(idx,1:3) accumarray(subs,mst{2}(:,4),[],@max)];
c = c+1;
save(fullfile(dd,sprintf('T%04d.mat',c)),'data','mst','-v7.3')
fprintf('%-40s%s\n',sprintf('T%04d.mat',c),datestr(now,'dd HH:MM:SS'))
        
% Set data to read only and hidden
fileattrib(fullfile(dd,'*.mat'),'+h -w','','s')

%% Generate grand master.mat file [replace with make master file]
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
[tck,~,posGrow] = unique(tck);
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
        
        idxM = mst{2}(:,1) == find(strcmpi(regexp(d(f).name,'\w+(?=_)','match'), tck)) & ...
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
%% Avoid splits [Do the check for splits only, should not find any]
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
%% Sort mst with From [Replace with sortMst]
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
end