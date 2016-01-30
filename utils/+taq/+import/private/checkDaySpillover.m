function prob = checkDaySpillover(path2mat)
% CHECKDAYSPILLOVER Checks if same-day data spans multiple files

% list files
mstfiles = dir(fullfile(path2mat,'*.mst'));
mstfiles = {mstfiles.name}';
nfiles   = numel(mstfiles);

% preallocate
tickers = table(cell(0,1),zeros(0,1,'uint32'),'VariableNames',{'Symb','Date'});
prob    = table(cell(0,1),zeros(0,1,'uint32'),cell(0,1), 'VariableNames',{'Symb','Date','File'});

h           = waitbar(0,'Checking...');
cleanupWait = onCleanup(@() delete(h));

map = containers.Map('keyType','char','ValueType','uint64');
for ii = 1:nfiles
    f = fullfile(path2mat,mstfiles{ii});
    s = load(f,'-mat');
    if isa(s.mst,'dataset'),
        s.mst = dataset2table(s.mst);
    end

    map = addKeyVal(map, s.ids);

    % Drop out-of-date
    ikeep   = min(s.mst.Date) - tickers.Date < 10;
    tickers = tickers(ikeep,:);

    % Check spillover
    s.mst.Symb = s.ids(s.mst.Id);
    keyA       = getValue(map,s.mst.Symb)*1e8 + uint64(s.mst.Date);
    keyB       = getValue(map,tickers.Symb)*1e8 + uint64(tickers.Date);
    [idx,pos]  = ismember(keyA,keyB);

    if any(idx)
        warning('Spillover: some data for the same day is spread over %s and %s.',mstfiles{ii-1}, mstfiles{ii})
        tmp                    = s.mst(idx,{'Symb','Date'});
        tmp.File               = repmat(mstfiles(ii),nnz(idx),1);
        prob                   = [prob; tmp];
        tickers.Date(pos(idx)) = s.mst.Date(idx);
    end

    % Add new
    tickers = [tickers; s.mst(~idx,{'Symb','Date'})];

    if mod(ii,100)
        waitbar(ii/nfiles,h)
    end
end
waitbar(1,h,'Completed')
pause(0.35)
end

function map = addKeyVal(map, keys)
inew = ~map.isKey(keys);
if any(inew)
    vals   = uint64(1:nnz(inew))+map.Count;
    map    = [map; containers.Map(keys(inew), vals)];
end
end

function val = getValue(map, keys)
if isempty(keys)
    val = zeros(0,1,'uint64');
else
    val = values(map,keys);
    val = cat(1,val{:});
end
end