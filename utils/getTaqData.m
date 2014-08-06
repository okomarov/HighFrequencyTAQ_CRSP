function out = getTaqData(master, tickers, from, to, varnames, path2data, debug)

% Checks and defaults
narginchk(2,7)
if isstring(tickers), tickers = {tickers}; end
if nargin < 3 || isempty(from), from = 0; end
if nargin < 4 || isempty(to),   to = inf; end
if nargin < 5 || isempty(varnames)
    varnames = '';
elseif isstring(varnames)
    varnames = {varnames};
end
if nargin < 6 || isempty(path2data),  path2data  = '.\data\TAQ'; end
if nargin < 7 || isempty(debug), debug = false; end   

    
% Find tickers in the master file
[ifound,ids] = ismember(upper(tickers),upper(master.ids));
if all(~ifound)
    warning('None of the TICKERS were found.')
    return
elseif any(~ifound)
    warning('\nThe following TICKERS were not matched:%s%s.',sprintf(' ''%s'',',tickers{~ifound}),char(8))
end

% Select entries matching the id and dates
master = master.mst(ismember(master.mst.Id,ids),:);
master = master(master.Date <= to & master.Date >= from,:);

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
    varnames   = unique(['Datetime' varnames],'stable');
else
    varnames   = unique(['Time' varnames],'stable');
end

% Display waitbar
h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
set(findall(h,'Type','text'),'interpreter','none')
setappdata(h,'canceling',false)

% Data files to load
files   = unique(master.File);
nfiles  = numel(files);
out     = cell(nfiles,1);
elapsed = zeros(nfiles+1,1);

% Open matlabpool
% poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'},'debug',debug)

tic
% LOOP all .mat files
% par
for ii = 1:nfiles
    % Check if progress was cancelled
    if getappdata(h,'canceling')
        break
    end
    
    % Update waitbar
    matname  = matnames{files(ii)};
    x        = (ii-1)/nfiles;
    time2end = elapsed(ii)*(1-x)/x;
    waitbar(x,h,sprintf('Loading file %s - remaining time %s',matname, sec2time(time2end)))
    
    % Load .mat file
    s = load(fullfile(path2data, matname),'data');
            
    % Records within the loaded data file
    mstfile = master(master.File == files(ii),:);
    idata   = mcolon(mstfile.From,mstfile.To);
    blocks  = mstfile.To - mstfile.From + 1;
    out{ii} = table(RunLength(mstfile.Id, blocks),'VariableNames',{'Id'}); 
    
    % Retrieve data
    if ~any(idatetime)
        out{ii}.Datetime = RunLength(yyyymmdd2serial(mstfile.Date),blocks) + hhmmssmat2serial(s.data.Time(idata,:));
    end
    for v = 1:numel(varnames)
        fname = varnames{v};
        out{ii}.(fname) = s.data.(fname)(idata,:);
    end

    % Update waitbar 
    waitbar(ii/nfiles,h)
    elapsed(ii+1) = toc;
end
out = cat(1,out{:});
delete(h)
end