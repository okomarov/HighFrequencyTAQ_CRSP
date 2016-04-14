function [out,ifound] = getTrades(idtype, id, date, path2data, updatebar)
% [out,ifound] = getTrades(idtype, id, date, path2data, updatebar)

% Checks and defaults
if nargin < 1,                          idtype    = 'symbol';       end
if nargin < 2,                          id        = [];             end
if nargin < 3 || isempty(date),         date      = inf;            end
if nargin < 4 || isempty(path2data),    path2data = '.\data\TAQ';   end
if nargin < 5,                          updatebar = false;          end
if isrowchar(id),  id  = {id}; end

[files_id, ifound, id] = filterFilesByID(idtype,id,path2data);

% Check if found matches
if all(~ifound)
    warning('None of the IDs were found.')
    out = [];
    return
elseif any(~ifound)
    warning('Some IDs were not matched. Check IFOUND output index.')
end

[files_date, date, dttype] = filterFilesByDate(date,path2data);

files  = intersect(files_id,files_date);
nfiles = numel(files);

% List files
dd       = dir(fullfile(path2data,'*.mat'));
matnames = {dd.name};
dd       = dir(fullfile(path2data,'*.idx'));
idxnames = {dd.name};

updatebar = updatebar && nfiles > 1;
out       = cell(nfiles,1);
elapsed   = 0;

% Display waitbar
if updatebar
    h              = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
    set(findall(h,'Type','text'),'interpreter','none')
    setappdata(h,'canceling',false)
    cleanupWaitbar = onCleanup(@()delete(h));
end

tic
for ii = 1:nfiles
    % Check if progress was cancelled
    if updatebar && getappdata(h,'canceling')
        break
    end

    mat_name = matnames{files(ii)};
    idx_name = idxnames{files(ii)};

    % Update waitbar
    if updatebar
        if ii == 1
            waitbar(0,h,sprintf('Loading file %s - ETA calculating',mat_name))
        else
            x        = (ii-1)/nfiles;
            time2end = elapsed*(nfiles-ii+1)/(ii-1); % (1-x)/x
            waitbar(x,h,sprintf('Loading file %s - ETA %s',mat_name, sec2time(time2end)))
        end
    end

    % Load data and index
    [index,symbol] = getIndexRecords(fullfile(path2data, idx_name), idtype, id, dttype, date);
    if isempty(index)
        continue
    end
    load(fullfile(path2data, mat_name));

    % Select data
    idata = mcolonint(index.From,index.To);
    data  = data(idata,:);

    % Add numeric id to symbol
    blocks = double(index.To - index.From + 1);
    if strcmpi(idtype,'Symbol')
        [~,pos]      = ismember(symbol(index.Id), id);
        data.Id(:,1) = uint32(RunLength(pos, blocks));
    end

    % Add date
    try
        data.Date(1);
    catch
        data.Date(:,1) = uint32(RunLength(index.Date, blocks));
    end

    out{ii} = data;

    % Update waitbar
    if updatebar
        waitbar(ii/nfiles,h)
        elapsed = toc;
    end
end
out = cat(1,out{:});
end

function [files, iskey, id] = filterFilesByID(idtype,id,path2data)
% Select the files that contain the queried IDs

if isempty(id)
    files = [];
    iskey = [];
else
    id2files = getRelevantIndexMap(idtype,path2data);
    if strcmpi(idtype,'symbol')
        id = upper(id);
    end
    iskey = id2files.isKey(id);
    files = values(id2files, id(iskey));
    files = unique([files{:}]);
end
end

function [files,date,dttype] = filterFilesByDate(date,path2data)
% Select the files that contain the queried DATEs

date2files = getRelevantIndexMap('date',path2data);
szDate     = size(date);
iinf       = isinf(date);

% [inf, inf] case
if isrow(date) && szDate(2) == 2 && all(iinf)
    date = inf;
end
% Inf case
if isscalar(date) && iinf
    dttype = 'all';
    files  = date2files.values;
    files  = unique([files{:}]);
    return
end

% Single date
if isscalar(date)
    dttype = 'scalar';
    dates  = {date};

% [from, to]
elseif szDate(2) == 2 && all(~iinf)
    dttype = 'fromto';
    dates  = num2cell(date(1):date(2));
    dates  = dates(date2files.isKey(dates));

% [inf, to] or [from, inf] cases
elseif szDate(2) == 2 && any(iinf)
    keys  = date2files.keys;
    keys  = [keys{:}];
    files = date2files.values;

    if iinf(1)
        dttype = 'le';
        files  = files(keys <= date(2));
    elseif iinf(2)
        dttype = 'ge';
        files  = files(keys >= date(1));
    end
    files = unique([files{:}]);

    return

% Column of dates
elseif szDate(2) == 1
    dttype = 'set';
    dates  = num2cell(date);
else
    error('DATE cannot have more than 2 columns.')
end

% Retrieve file list
files = values(date2files, dates);
files = unique([files{:}]);
end

function map = getRelevantIndexMap(type,path2data)
% Load the maps that associate ID/DATE with list of files

PREFIX_INDEX = 'index_';

if any(strcmpi(type, {'symbol','id','permno','date'}))
    indexname = fullfile(path2data, sprintf([PREFIX_INDEX '%s'],type));
    s         = load(indexname,'-mat');
    map       = s.(type);
else
    error('getTaqData:invalidIdtype','IDTYPE can be ''symbol'', ''permno'' or ''id''.')
end
end

function [index,symbol] = getIndexRecords(filename, idtype, id, dttype, date)
% Get master records that meet ID and DATE criteria

s      = load(filename,'-mat');
index  = s.index;
symbol = s.symbol;

switch dttype
    case 'scalar'
        idx = index.Date == date;
    case 'all'
        idx = true(size(index,1),1);
    case 'ge'
        idx = index.Date >= date(1);
    case 'le'
        idx = index.Date <= date(2);
    case 'fromto'
        idx = index.Date >= date(1) & index.Date <= date(2);
    case 'set'
        idx = ismember(index.Date, date);
end
index = index(idx,:);

switch idtype
    case 'symbol'
        id  = find(ismember(symbol,id));
        idx = ismember(index.Id, id);
end
index = index(idx,:);
end
