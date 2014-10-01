function out = getTaqData(master, tickers, from, to, varnames, path2data)

% GETTAQDATA Retrieve TAQ data from a given directory
% 
%   GETTAQDATA(MASTER, TICKERS, FROM, TO, VARNAMES, PATH2DATA)
%
% 
% Example:
%   
%   path2data = '.\data\TAQ\sampled';
%   master    = load(fullfile(path2data, 'master'), '-mat');
%   data      = getTaqData(master, {'A','AA'}, 19961231, 20080122, 'Price', path2data);

addpath .\utils\mcolon

% Checks and defaults
narginchk(1,6)
if nargin < 2, tickers = []; end
if isstring(tickers), tickers = {tickers}; end
if nargin < 3 || isempty(from), from = 0; end
if nargin < 4 || isempty(to),   to = inf; end
if nargin < 5 || isempty(varnames)
    varnames = '';
elseif isstring(varnames)
    varnames = {varnames};
end
if nargin < 6 || isempty(path2data),  path2data  = '.\data\TAQ'; end
% if nargin < 7 || isempty(debug), debug = false; end   

    
% Find tickers in the master file
if ~isempty(tickers)
    [ifound,ids] = ismember(upper(tickers),upper(master.ids));
    if all(~ifound)
        warning('None of the TICKERS were found.')
        return
    elseif any(~ifound)
        warning('\nThe following TICKERS were not matched:%s%s.',sprintf(' ''%s'',',tickers{~ifound}),char(8))
    end
    
    % Filter based on IDs
    master = master.mst(ismember(master.mst.Id,ids),:);
end

% Filter based on dates
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

% Type of IDs available
hasUnID   = any(strcmpi(getVariableNames(master), 'UnID'));
hasPermno = any(strcmpi(getVariableNames(master), 'permno'));

% Data files to load
files     = unique(master.File);
nfiles    = numel(files);
updatebar = nfiles > 1;
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
    matname  = matnames{files(ii)};
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
        out{ii}.Datetime = RunLength(yyyymmdd2serial(mstfile.Date),blocks) + hhmmssmat2serial(s.data.Time(idata,:));
    end
    for v = 1:numel(varnames)
        fname = varnames{v};
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