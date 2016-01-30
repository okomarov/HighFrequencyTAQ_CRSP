function [out,ifound] = getTrades(idtype, id, date, path2data, updatebar)


% Checks and defaults
if nargin < 1,                          idtype    = 'symbol';       end
if nargin < 2,                          id        = [];             end
if nargin < 3 || isempty(date),         date      = inf;            end
if nargin < 4 || isempty(path2data),    path2data = '.\data\TAQ';   end
if nargin < 5,                          updatebar = true;           end
if isrowchar(id),  id  = {id}; end

% Filter by ID
if isempty(id)
    filesId = [];
else
    
    switch lower(idtype)
        case 'symbol'
            % Load master
            load(fullfile(path2data,'master_symbol'),'-mat');
            filesId = values(mstSymb, upper(id));
            filesId = unique([filesId{:}]);

        case 'permno'
            % Load master
            load(fullfile(path2data,'master_permno'),'-mat');
            filesId = values(mstPermno, id);
            filesId = unique([filesId{:}]);
            
        case 'id'
            % Load master
            load(fullfile(path2data,'master_symbol'),'-mat');
            filesId = mstSymb.File(id);
            filesId = unique([filesId{:}]);

        otherwise
            error('getTaqData:invalidIdtype','IDTYPE can be ''symbol'', ''permno'' or ''id''.')
    end

    % Check if found matches
    if isempty(filesId)
        warning('None of the IDs were found.')
        out = [];
        return
    elseif any(~ifound)
        warning('The following IDs were not matched:%s%s.',sprintf(' ''%s'',',id{~ifound}),char(8))
    end
end

load(fullfile(path2data,'master_date'),'-mat');
szDate = size(date);
if isrow(date) && szDate(2) == 2 && all(isinf(date))
    date = inf;
end
if isscalar(date) && isinf(date)
    filesDate = mstDate.values;
else
    if szDate(2) == 1
        % Do nothing here
    elseif szDate(2) == 2
        iinf = isinf(date);
        keys = mstDate.keys;
        keys = [keys{:}];
        if iinf(1)
            date = keys(keys <= date(2));
        elseif iinf(2)
            date = keys(keys >= date(1));
        else
            dates = date(1):date(2);
            date = dates(mstDate.isKey(dates));
        end
    else
        error('DATE cannot have more than 2 columns.')
    end
end
filesDate = values(mstDate,num2cell(date));

% Filter based on dates
if from ~= 0 || isfinite(to)
    master = master(in(master.Date,[from, to]),:);
end
end