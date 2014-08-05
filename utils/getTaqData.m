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
    buildtable = '@(r) table(zeros(r,1,''uint16'')';
else
    varnames   = unique(['Time' varnames],'stable');
    buildtable = '@(r) table(zeros(r,1,''uint16''), NaN(r,1)';
end

% Create pre-allocation function for table construction
for v = 1:numel(varnames)
    field = tmp.data.(varnames{v});
    ncols = size(field,2);
    if isinteger(field)
        buildtable = [buildtable sprintf(', zeros(r,%d,''%s'')',ncols,class(field))];
    elseif isfloat(field)
        buildtable = [buildtable sprintf(', NaN(r,%d,''%s'')',ncols,class(field))];
    else
        buildtable = [buildtable sprintf(', repmat('' '',r, %d)',ncols)];
    end
end
buildtable = [buildtable, ', ''VariableNames'', {''Id''', sprintf(',''%s''', varnames{:}), '})'];
buildtable = str2func(buildtable);

% % Display waitbar
% h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
% set(findall(h,'Type','text'),'interpreter','none')
% setappdata(h,'canceling',false)

% Data files to load
files   = unique(master.File);
nfiles  = numel(files);
res     = cell(nfiles,1);
% elapsed = zeros(nfiles+1,1);

% Open matlabpool
% poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'},'debug',debug)

% tic
% LOOP all .mat files
% par
for ii = 1:nfiles
%     % Check if progress was cancelled
%     if getappdata(h,'canceling')
%         break
%     end
    
    % Update waitbar
    matname  = matnames{files(ii)};
%     x        = (ii-1)/nfiles;
%     time2end = elapsed(ii)*(1-x)/x;
%     waitbar(x,h,sprintf('Loading file %s - remaining time %s',matname, sec2time(time2end)))
    
    % Load .mat file
    s = load(fullfile(path2data, matname),'data');
            
    % Records within the loaded data file
    mstfile       = master(master.File == files(ii),:);
    [mstsqz,isrt] = squeezeFromTo(mstfile(:,{'Id','From','To'}));
    dates         = yyyymmdd2serial(mstfile.Date(isrt));
    nrec          = size(mstsqz,1);
    
    % Count how many datapoints per master record
    blocks  = mstsqz.To - mstsqz.From + 1;
    cblocks = [0; cumsum(blocks)];
    
    % Preallocate
    res{ii} = buildtable(cblocks(end));

    % Assign Id
    res{ii}.Id = RunLength(mstsqz.Id, blocks);
    
    % Assign date
    res{ii}.Datetime = RunLength(dates,  mstfile.To - mstfile.From + 1);
    
    % Loop for each record
    for jj = 1:nrec
        % Create indices pointing to the output and input
        iout  = cblocks(jj)+1:cblocks(jj+1);
        idata = mstsqz.From(jj):mstsqz.To(jj);

        % Retrieve data
        if ~any(idatetime)
            res{ii}.Datetime(iout) = res{ii}.Datetime(iout) + hhmmssmat2serial(s.data.Time(idata,:));
        end
        for v = 1:numel(varnames)
            fname = varnames{v};
            res{ii}.(fname)(iout,:) = s.data.(fname)(idata,:);
        end
    end
%     % Update waitbar 
%     waitbar(ii/nfiles,h)
%     elapsed(ii+1) = toc;
end


% if getappdata(h,'canceling')
%     out = out(1:find(out.Id == 0,1,'first')-1,:);
% end
% delete(h)
end