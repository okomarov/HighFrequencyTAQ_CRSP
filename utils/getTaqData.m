function [out,ifound] = getTaqData(idtype, id, from, to, varnames, path2data, master, updatebar)

% GETTAQDATA Retrieve TAQ data from a given directory
% 
%   GETTAQDATA(IDTYPE, ID, FROM, TO, VARNAMES, PATH2DATA)
%
% 
% Example:
%   
%   path2data = '.\data\TAQ\sampled\5min';
%   data      = getTaqData('symbol', {'A','AA'}, 19961231, 20080122, 'Price', path2data);

% Checks and defaults
if nargin < 1,                          idtype    = 'symbol';       end
if nargin < 2,                          id        = [];             end
if nargin < 3 || isempty(from),         from      = 0;              end
if nargin < 4 || isempty(to),           to        = inf;            end
if nargin < 5,                          varnames  = [];             end
if nargin < 6 || isempty(path2data),    path2data = '.\data\TAQ';   end
if nargin < 7 || isempty(master)
    try
        master = load(fullfile(path2data, 'master'), '-mat');
    catch
        master    = load(fullfile('..',path2data, 'master'), '-mat');
        path2data = fullfile('..',path2data);
    end
end
if nargin < 8,                          updatebar = true;           end
if isstring(id),  id  = {id}; end
if isstring(varnames), varnames = {varnames}; end

% Filter by ID
if isempty(id)
    try
        master = master.mst;
    catch
    end
else
    
    % Find tickers in the master file
    switch lower(idtype)
        case 'symbol'
            [ifound, ids] = ismember(upper(id), upper(master.ids));
            imst          = ismember(master.mst.Id,ids);
            
        case 'permno'
            try 
                id2permno = loadresults('masterPermno');
            catch 
                id2permno = loadresults('masterPermno','..\results');
            end
            idx       = ismember(id2permno.Permno,id);
            id2permno = id2permno(idx,:);
            ifound    = ismember(id, id2permno.Permno);
            imst      = ismembIdDate(master.mst.Id, master.mst.Date,id2permno.Id, id2permno.Date);
            id        = cellstr(num2str(id(:)));
                    
        case 'id'
            imst   = ismember(master.mst.Id, id);
            ifound = ismember(id, master.mst.Id(imst));
            id     = cellstr(num2str(id(:)));
            
        otherwise
            error('getTaqData:invalidIdtype','IDTYPE can be ''symbol'', ''permno'' or ''id''.')
    end

    % Check if found matches
    if all(~ifound)
        warning('None of the IDs were found.')
        out = [];
        return
    elseif any(~ifound)
        warning('The following IDs were not matched:%s%s.',sprintf(' ''%s'',',id{~ifound}),char(8))
    end
    
    % Filter based on IDs
    master = master.mst(imst,:);
end

% Filter based on dates
if from ~= 0 || isfinite(to)
    master = master(in(master.Date,[from, to]),:);
end

% List files
dd       = dir(fullfile(path2data,'*.mat'));
matnames = {dd.name};

% Retrieve varnames
tmp            = load(fullfile(path2data, matnames{end}),'data');
availablenames = getVariableNames(tmp.data);
if isempty(varnames)
    varnames = setdiff(availablenames,'Time','stable');
end

% Check if Datetime already available
idatetime = any(strcmpi(availablenames, 'Datetime'));
if idatetime
    varnames = unique(['Datetime' varnames],'stable');
else
    varnames = unique(['Time' varnames],'stable');
end

% Type of IDs available
hasUnID   = any(strcmpi(getVariableNames(master), 'UnID'));
hasPermno = any(strcmpi(getVariableNames(master), 'permno'));

% Data files to load
files     = unique(master.File);
nfiles    = numel(files);
updatebar = updatebar && nfiles > 1;
out       = cell(nfiles,1);
elapsed   = zeros(nfiles+1,1);

% Display waitbar
if updatebar
    h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
    set(findall(h,'Type','text'),'interpreter','none')
    setappdata(h,'canceling',false)
end
    
% Open matlabpool
% poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'},'debug',debug)

tic
% LOOP all .mat files
% par
for ii = 1:nfiles
    % Check if progress was cancelled
    if updatebar && getappdata(h,'canceling')
        break
    end
    
    % Update waitbar
    matname = matnames{files(ii)};
    if updatebar
        x        = (ii-1)/nfiles;
        time2end = elapsed(ii)*(1-x)/x;
        waitbar(x,h,sprintf('Loading file %s - remaining time %s',matname, sec2time(time2end)))
    end
    
    % Load .mat file
    s = load(fullfile(path2data, matname),'data');
            
    % Records within the loaded data file
    mstfile = master(master.File == files(ii),:);
    idata   = mcolon(mstfile.From,mstfile.To);
    blocks  = double(mstfile.To - mstfile.From + 1);
    Id      = RunLength(mstfile.Id, blocks);
    out{ii} = table(Id(:),'VariableNames',{'Id'}); 
    if hasUnID
        out{ii}.UnID = reshape(RunLength(mstfile.UnID, blocks),[],1);
    end
    if hasPermno
        out{ii}.Permno = reshape(RunLength(mstfile.Permno, blocks),[],1);
    end
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
if updatebar
    delete(h)
end
end

function idx = mcolon(from, to)
nobs     = to-from+1;
idx      = ones(sum(nobs),1);
pos      = cumsum(nobs(1:end-1))+1;
df       = from(2:end) - to(1:end-1);
idx(pos) = df;
idx(1)   = from(1);
idx      = cumsum(idx);
end