function [out,ifound] = getTrades(idtype, id, date, path2data, updatebar)


% Checks and defaults
if nargin < 1,                          idtype    = 'symbol';       end
if nargin < 2,                          id        = [];             end
if nargin < 3 || isempty(date),         date      = inf;            end
if nargin < 4 || isempty(path2data),    path2data = '.\data\TAQ';   end
if nargin < 5,                          updatebar = true;           end
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

files     = intersect(files_id,files_date);
nfiles    = numel(files);

% List files
dd       = dir(fullfile(path2data,'*.mat'));
matnames = {dd.name};
dd       = dir(fullfile(path2data,'*.mst'));
mstnames = {dd.name};

updatebar = updatebar && nfiles > 1;
out       = cell(nfiles,1);
elapsed   = zeros(nfiles+1,1);

% Display waitbar
if updatebar
    h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
    set(findall(h,'Type','text'),'interpreter','none')
    setappdata(h,'canceling',false)
    cleanupWaitbar = onCleanup(@()delete(h));
end

for ii = 1:nfiles
    % Check if progress was cancelled
    if updatebar && getappdata(h,'canceling')
        break
    end

    % Update waitbar
    matname = matnames{files(ii)};
    if updatebar
        if ii == 1
            waitbar(0,h,sprintf('Loading file %s - ETA calculating',matname))
        else
            x        = (ii-1)/nfiles;
            time2end = elapsed(ii)*(nfiles-ii+1)/(ii-1); % (1-x)/x
            waitbar(x,h,sprintf('Loading file %s - ETA %s',matname, sec2time(time2end)))
        end
    end

    % Load .mat file
    s = load(fullfile(path2data, matname));

    % Records within the loaded data file
    mstfile = master(master.File == files(ii),:);
    idata   = mcolonint(mstfile.From,mstfile.To);
    blocks  = double(mstfile.To - mstfile.From + 1);
    Id      = RunLength(mstfile.Id, blocks);
    out{ii} = table(Id(:),'VariableNames',{'Id'});
    
    % Retrieve data
    if ~any(idatetime)
        dates            = RunLength(yyyymmdd2serial(mstfile.Date),blocks);
        out{ii}.Datetime = dates(:) + hhmmssmat2serial(s.data.Time(idata,:));
    end
    for v = 1:numel(varnames)
        fname           = varnames{v};
        out{ii}.(fname) = s.data.(fname)(idata,:);
    end

    % Update waitbar
    if updatebar
        waitbar(ii/nfiles,h)
        elapsed(ii+1) = toc;
    end
end
out = cat(1,out{:});
end

function [files, iskey, id] = filterFilesByID(idtype,id,path2data)
if isempty(id)
    files = [];
    iskey = [];
else
    id2files = getRelevantMasterMap(idtype,path2data);
    if strcmpi(idtype,'symbol')
        id = upper(id);
    end
    iskey = id2files.isKey(id);
    files = values(id2files, id(iskey));
    files = unique([files{:}]);
end
end

function [files,dates,dttype] = filterFilesByDate(date,path2data)
date2files = getRelevantMasterMap('date',path2data);
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
    dates  = inf;
    return
end

% Single date
if isscalar(date)
    dttype = 'scalar';
    dates  = date;
% [from, to]     
elseif szDate(2) == 2 && all(~inf)
    dttype = 'set';
    dates  = date(1):date(2);
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
    dates = date;
    return
% Column of dates
elseif szDate(2) == 1
    dttype = 'set';
    dates  = date;
else
    error('DATE cannot have more than 2 columns.')
end

% Retrieve file list
files = values(date2files, num2cell(dates));
files = unique([files{:}]);
end

function map = getRelevantMasterMap(type,path2data)
if any(strcmpi(type, {'symbol','id','permno','date'}))
    mstname = fullfile(path2data, sprintf('master_%s',type));
    s       = load(mstname,'-mat');
    map     = s.(type);
else
    error('getTaqData:invalidIdtype','IDTYPE can be ''symbol'', ''permno'' or ''id''.')
end
end

