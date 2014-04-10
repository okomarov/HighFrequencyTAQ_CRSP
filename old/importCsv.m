%% FINAL files format

% 'data', contains trade data in a 1x6 cell:
%
% 1     [hh, mm, ss]                            uint16
% 2     price                                   single
% 3     size                                    uint32
% 4     [G127 rule,  correction]                uint16
% 5     condition                               2char
% 6     exchange                                1char

% 'mst', master file with ticker symbol, dates, and starting positions:
%
% 1     tickers                                 cellstr
% 2     [idxTickers, yyyymmdd, startPosData]	uint32
% 3     startPosMst                             double      length as tickers+1 because has end+1 position

%% Import CSV
dd   = 'E:\data';
%d    = struct('name','b045446dd7078060.csv');
d    = dir(fullfile(dd,'raw','*.zip'));
opt  = {'Delimiter',',','HeaderLines',1};
N    = 1e5;
name = '';
% Loop each csv
for f = 1:numel(d)
    filename = unzip(fullfile(dd,'raw',d(f).name),dd);
    %filename = {'E:\data\b045446dd7078060.csv'};
    fid      = fopen(filename{:});
    fprintf('%-40s%s\n',filename{:},datestr(now,'dd HH:MM:SS'))
    
    % Reset filepart counter if previous filename was different
    if ~strcmp(d(f).name(1:6), name)
        c = 0;
    end
    
    while ~feof(fid)
        c    = c + 1;
        % Approximate preallocation for 1.6GB
        nite = 200;
        data = cell(nite,11);
        mst  = cell(nite,2);
        ii   = 0;
        while ii < nite && ~feof(fid)
            ii = ii + 1;
            %    1       2       3-5       6       7        8           9           10          11
            % ticker | date | hh:mm:ss | price | size | G127 rule | correction | condition | exchange
            data(ii,:)         = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',N,opt{:});
            % Keep copy for check (redundant but faster)
            cmp                = data(ii,:);
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
        end
        
        % Collect everything in one cell per column data and regorganize a bit
        data            = arrayfun(@(x) cat(1,data{:,x}),3:11,'un',0);
        data{1}         = cat(2,data{1:3});
        data{6}         = cat(2,data{6:7});
        data([2:3 7])   = [];
        mst             = {cat(1,mst{:,1}) cat(1,mst{:,2})}; % cellfun('length',mst(1:ii-1,1))};
        % Translate block indexing to whole file indexing
        idx             = [false; mst{2}(2:end,3) == 1];     % Start position of each block is 1 (trivially)
        val             = N-mst{2}(find(idx)-1,3)+1;         % Calculate within blocks last-position increments
        mst{2}(:,[1,3]) = [1 1; diff(mst{2}(:,[1, 3]))];     % Diff positions to get increments
        mst{2}(idx,1)   = 1;                                 % Reset 0s (recall it's uint type) to 1
        mst{2}(idx,3)   = val;                               % Set all last increments within block
        mst{2}(:,[1 3]) = cumsum(mst{2}(:,[1,3]));           % Cumulate back to get file positions
        % Remove ticker repetitions and redundant ticker-date pairs retaining only first start pos
        [mst{1},~,pos]  = unique(mst{1});
        mst{2}(:,1)     = pos(mst{2}(:,1));
        [~,idx]         = unique(mst{2}(:,1:2),'rows','first');
        mst{2}          = mst{2}(idx,:);
        % Start positions for each ticker within mst
        mst{3}          = find([1; diff(mst{2}(:,1)); 1]);
        %fprintf('%-40s%s\n','formatted',datestr(now,'dd HH:MM:SS'))
        
        % Saving cell data and master as is - tested 1.6GB - 21 sec to save; 4.5 sec to load; -30% space
        name = d(f).name(1:6); % save it for comparison later
        save(fullfile(dd,sprintf('%s_%02d.mat',name,c)),'data','mst','-v7.3')
        fprintf('%-40s%s\n',sprintf('%s_%02d.mat',name,c),datestr(now,'dd HH:MM:SS'))
        
        % Chesck .mat against csv
        len      = numel(cmp{1});
        [~,~,id] = unique(cmp{1});
        [un,pos] = unique([id cmp{2}],'rows','first');
        tests    = [...
        isequal(diff([un(2:end,:) pos(2:end)+ N*(ii-1)],[],1), diff(mst{2}(mst{2}(:,3) > pos(1)*N*(ii-1),:),[],1)) || numel(pos)==1
        isequal(cat(2,cmp{3:5}), data{1}(end-len+1:end,:))
        isequal(         cmp{6}, data{2}(end-len+1:end,:))
        isequal(         cmp{7}, data{3}(end-len+1:end,:))
        isequal(cat(2,cmp{8:9}), data{4}(end-len+1:end,:))
        isequal(  char(cmp{10}), data{5}(end-len+1:end,:)) || isequal(char(cmp{10}), data{5}(end-len+1:end,1))
        isequal(  char(cmp{11}), data{6}(end-len+1:end,:))];
        fprintf('%-40s%s\n',sprintf('OK: %d/7',sum(tests)),datestr(now,'dd HH:MM:SS'))
    end
    fclose(fid);
    pause(0.5)
    delete(filename{:})
end

% Set data to read only and hidden
fileattrib(fullfile(d,'*.mat'),'+h -w','','s')