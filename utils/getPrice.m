function out = getPrice(master, tickers, from, to)

narginchk(2,4)
if nargin < 3 || isempty(from), from = 0; end
if nargin < 4 || isempty(to),   to = inf; end
if ischar(tickers) && isrow(tickers)
    tickers = {tickers};
end
    
% Find tickers in the master file
[ifound,ids] = ismember(tickers,master.ids);
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

% Preallocate output
nrows = cblocks(end);
out   = table(zeros(nrows,1,'uint16'),... 
               NaN(nrows,1),...
               NaN(nrows,1,'single'),...
               'VariableNames',{'Id','Datetime','Price'});

% Display waitbar
h = waitbar(0,'','CreateCancelBtn','setappdata(gcbf,''canceling'',true)');
setappdata(h,'canceling',false)

% Data files to load
files   = unique(master.File);
nfiles  = numel(files);

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
    matname = sprintf('T%04d.mat',files(ii));
    x = (ii-1)/nfiles;
    time2end = elapsed(ii)*(1-x)/x;
    waitbar(x,h,sprintf('Loading file %s - remaining time %s',matname, sec2time(time2end)))
    
    % Load .mat file
    s = load(fullfile('.\data\TAQ', matname),'data');
            
    % Records within the loaded data file
    mstfile = master(master.File == files(ii),:);
    dates   = yyyymmdd2serial(double(mstfile.Date));
    nrec    = size(mstfile,1);
    
    % Loop for each record
    for jj = 1:nrec
        % Create indices pointing to the output and input
        c     = c+1;
        iout  = cblocks(c)+1:cblocks(c+1);
        idata = mstfile.From(jj):mstfile.To(jj);
        % Retrieve data
        out.Id      (iout) = mstfile.Id(jj);
        out.Datetime(iout) = dates(jj) + hhmmssmat2serial(s.data.Time(idata,:));
        out.Price   (iout) = s.data.Price(idata);
        
        % Update waitbar 
        waitbar(x + jj/(nfiles*nrec),h)
    end
    elapsed(ii+1) = toc;
end
if getappdata(h,'canceling')
    out = out(1:find(out.Id == 0,1,'first')-1,:);
end
delete(h)
end