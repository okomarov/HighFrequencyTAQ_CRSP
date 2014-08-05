function out = getData(master, tickers, from, to, varnames, path2data)
narginchk(2,6)

if nargin < 3 || isempty(from), from = 0; end
if nargin < 4 || isempty(to),   to = inf; end
if nargin < 5 || isempty(varnames)
    varnames = '';
elseif isstring(varnames)
    varnames = {varnames};
end
if nargin < 6 || isempty(path2data),  path2data  = '.\data\TAQ'; end
   
if isstring(tickers)
    tickers = {tickers};
end
    
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

% Count how many datapoints per master record
blocks  = master.To - master.From + 1;
cblocks = [0; cumsum(blocks)];

% List files
dd       = dir(fullfile(path2data,'*.mat'));
matnames = {dd.name};

% Retrieve varnames
tmp    = load(fullfile(path2data, matnames{end}),'data');
vnames = getVariableNames(tmp.data);
if isempty(varnames)
    varnames = setdiff(vnames,{'Time','Properties'},'stable');
end

% Check if it already has Datetime
idatetime = any(strcmpi(vnames, 'Datetime'));
if idatetime
    varnames = unique(['Datetime' varnames],'stable');
else
    varnames = unique(['Time' varnames],'stable');
end

% % Preallocate output
% nrows = cblocks(end);
% out   = table(zeros(nrows,1,'uint16'),... 
%                NaN(nrows,1),...
%                'VariableNames',{'Id','Datetime'});
% for ii = 1:numel(varnames)
%     field = varnames{ii};
%     if isinteger(tmp.data.(field))
%         out.(field) = zeros(nrows,size(tmp.data.(field),2),'like', tmp.data.(field));
%     elseif isfloat(tmp.data.(field))
%         out.(field) = NaN(nrows,size(tmp.data.(field),2),'like', tmp.data.(field));
%     else
%         out.(field) = repmat(' ', nrows,size(tmp.data.(field),2));
%     end
% end

% Display waitbar
h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
set(findall(h,'Type','text'),'interpreter','none')
setappdata(h,'canceling',false)

% Data files to load
files   = unique(master.File);
nfiles  = numel(files);

res     = cell(nfiles,1);

c = 0;
elapsed = zeros(nfiles+1,1);
tic
% LOOP all .mat files
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
    mstfile       = master(master.File == files(ii),:);
    [mstsqz,isrt] = squeezeFromTo(mstfile(:,{'From','To'}));
    dates         = yyyymmdd2serial(mstfile.Date(isrt));
    nrec          = size(mstsqz,1);
    
    % Preallocate
    res{ii} = table(zeros(nrec,1,'uint16'),... 
                    NaN(nrec,1),...
                    'VariableNames',{'Id','Datetime'});
    for v = 1:numel(varnames)
        fname = varnames{v};
        field = tmp.data.(fname);
        ncols = size(field,2);
        if isinteger(field)
            res{ii}.(fname) = zeros(nrec,ncols,'like', field);
        elseif isfloat(field)
            res{ii}.(fname) = NaN(nrec,ncols,'like', field);
        else
            res{ii}.(fname) = repmat(' ', ncols);
        end
    end
    
    % Loop for each record
    for jj = 1:nrec
        % Create indices pointing to the output and input
        c     = c+1;
        iout  = cblocks(c)+1:cblocks(c+1);
        idata = mstfile.From(jj):mstfile.To(jj);

        % Retrieve data
        out.Id(iout) = mstfile.Id(jj);
        if ~any(idatetime)
            out.Datetime(iout) = dates(jj) + hhmmssmat2serial(s.data.Time(idata,:));
        end
        for v = 1:numel(varnames)
            field = varnames{v};
            out.(field)(iout,:) = s.data.(field)(idata,:);
        end
    end
    % Update waitbar 
    waitbar(ii/nfiles,h)
    elapsed(ii+1) = toc;
end
if getappdata(h,'canceling')
    out = out(1:find(out.Id == 0,1,'first')-1,:);
end
delete(h)
end