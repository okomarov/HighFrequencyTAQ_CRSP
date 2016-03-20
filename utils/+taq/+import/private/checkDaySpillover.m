function prob = checkDaySpillover(path2mat)
% CHECKDAYSPILLOVER Checks if same-day data spans multiple files

% list files
indexfiles = dir(fullfile(path2mat,'*.idx'));
indexfiles = {indexfiles.name}';
nfiles     = numel(indexfiles);

% preallocate
tickers = table(cell(0,1),zeros(0,1,'uint32'),'VariableNames',{'Symb','Date'});
prob    = table(cell(0,1),zeros(0,1,'uint32'),cell(0,1), 'VariableNames',{'Symb','Date','File'});

h           = waitbar(0,'Checking...');
cleanupWait = onCleanup(@() delete(h));

map = containers.Map('keyType','char','ValueType','uint64');
for ii = 1:nfiles
    f = fullfile(path2mat,indexfiles{ii});
    s = load(f,'-mat');
    if isa(s.index,'dataset'),
        s.index = dataset2table(s.index);
    end

    map = addKeyVal(map, s.symbol);

    % Drop out-of-date
    ikeep   = min(s.index.Date) - tickers.Date < 10;
    tickers = tickers(ikeep,:);

    % Check spillover
    s.index.Symb = s.symbol(s.index.Id);
    keyA         = getValue(map,s.index.Symb)*1e8 + uint64(s.index.Date);
    keyB         = getValue(map,tickers.Symb)*1e8 + uint64(tickers.Date);
    [idx,pos]    = ismember(keyA,keyB);

    if any(idx)
        warning('Spillover: some data for the same day is spread over %s and %s.',indexfiles{ii-1}, indexfiles{ii})
        tmp                    = s.index(idx,{'Symb','Date'});
        tmp.File               = repmat(indexfiles(ii),nnz(idx),1);
        prob                   = [prob; tmp];
        tickers.Date(pos(idx)) = s.index.Date(idx);
    end

    % Add new
    tickers = [tickers; s.index(~idx,{'Symb','Date'})];

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
    vals = uint64(1:nnz(inew))+map.Count;
    map  = [map; containers.Map(keys(inew), vals)];
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