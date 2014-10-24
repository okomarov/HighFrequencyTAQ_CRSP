function [res, filename] = Analyze(fun, varnames, cached, path2mat, debug, varargin)

% ANALYZE Executes specified fun in parallel on the whole database (all .mat files)
%
%   ANALYZE(FUN, VARNAMES) FUN should a string with the name of one of
%                          the following sub-functions:
%                               - 'dailystats'
%                               - 'badprices'
%                               - 'avgtimestep'
%                          VARNAMES is a cell-array of strings (or string)
%                          with the VarNames of the dataset with the results
%
%   ANALYZE(..., PATH2DATA) If you wanna use other than '.\data\TAQ\T*.mat'
%                           files (default), then specify a different
%                           PATH2DATA with the initial pattern of the name,
%                           e.g '.\data\TAQ\sampled\S5m_*.mat'
%
%   ANALYZE(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the
%                        needed CACHED results.
%   ANALYZE(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be
%                       able to step through the code in debug mode.

if nargin < 2 || isempty(varnames);  varnames  = {'data','mst','ids'};  end
if nargin < 3,                       cached    = [];                    end
if nargin < 4 || isempty(path2mat);  path2mat  = '.\data\TAQ\T*.mat';    end
if nargin < 5 || isempty(debug);     debug     = false;                 end

% Simply call the specific subroutine
addpath(genpath('.\utils'))
writeto = '.\results\';
root    = fileparts(path2mat);

% Open matlabpool
poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'},'debug',debug)

% Get email credentials if not in debug
if ~debug; setupemail; end

fhandle = str2func(fun);

try
    tic
    dd  = dir(path2mat);
    N   = numel(dd);
    res = deal(cell(N,1));
    
    % Slice cached
    if isempty(cached)
        cached = cell(N,1);
    elseif ~iscell(cached)
        vnames = setdiff(getVariableNames(cached),'File','stable');
        cached = accumarray(cached.File,(1:size(cached))',[],@(x) {cached(x,vnames)});
    end
    
    % LOOP in parallel
    parfor f = 1:N
        disp(f)
        % Load data
        s      = load(fullfile(root,dd(f).name), varnames{:});
        cache  = [cached(f,:), f];
        % Convert to tables
        if isfield(s,'data') && isa(s.data,'dataset'), s.data = dataset2table(s.data); end
        if isfield(s,'mst')  && isa(s.mst ,'dataset'), s.mst  = dataset2table(s.mst ); end
        % Apply function
        res{f} = fhandle(s, cache, varargin{:});
    end
    % Collect all results and convert to dataset
    res = cat(1,res{:});
    
    % Export results and notify
    filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),fun);
    save(fullfile(writeto,filename), 'res')
    message = sprintf('Task ''%s'' terminated in %s',fun,sec2time(toc));
    disp(message)
    if ~debug, sendmail('o.komarov11@imperial.ac.uk', message,''); end
catch err
    filename = fullfile(writeto, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    if ~debug
        sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',fun), err.message, {filename})
    end
    rethrow(err)
end
if ~debug
    delete(gcp('nocreate'))
%     rmpref('Internet','SMTP_Password')
end
end
%% Subfunctions 

% Maximum trades per second
function res = maxtradepsec(s,cached)
% Count per id, date and second
Id          = cumsum([ones(1,'uint32'); diff(int8(s.data.Time(:,1))) < 0]);
Dates       = RunLength(s.mst.Date, double(s.mst.To-s.mst.From+1));
[un,~,subs] = unique([Id, Dates, s.data.Time],'rows');
counts      = [un accumarray(subs,1)];
% Pick max per date 
[date,~,subs] = unique(counts(:,2));
res           = table(date, accumarray(subs,counts(:,end),[],@max),'VariableNames',{'Date','Maxpsec'});
end

% Median price (net of the first selection step)
function res = medianprice(s,cached)
cachedmst = cached{1};

% STEP 1) Selection
inan = selecttrades(s.data);

% Prepare for daily median
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(cachedmst,1);
subs = uint32(RunLength((1:nmst)',nobs));

% Daily median price
cachedmst.MedPrice = accumarray(subs(~inan), s.data.Price(~inan),[nmst,1], @fast_median);
cachedmst.MedPrice(cachedmst.MedPrice == 0) = NaN;

res = cachedmst;
end

% Identify bad prices
function res = badprices(s, cached, dailycut)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices are < than .65x daily median or > than 1.5x daily median
[~,igoodprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);

% STEP 3) Bad days
res          = cached(:,{'Id','Date'});
subs         = uint32(RunLength((1:size(cached,1))',nobs));
res.Nbadsel  = uint32(accumarray(subs,  inan));
res.Nbadtot  = uint32(accumarray(subs,  inan | igoodprice ~= 1));
res.Isbadday = res.Nbadtot > ceil(dailycut*nobs);
end

% Count how many observations we loose in consolidation step
function res = consolidationcounts(s,cached)
cached = cached{1};

% Number of observations per day
nobs   = double(s.mst.To - s.mst.From + 1);
nmst   = size(s.mst,1);

% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,iprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
iprice     = iprice ~= 1;
inan       = inan | iprice;

% STEP 3) Bad series/days
inan = inan | RunLength(cached.Isbadday,nobs);

% STEP 4) Count how many observations we loose from median consolidation
mstrow     = RunLength((1:nmst)',nobs);
% Count by mst row and timestamp
[un,~,subs] = unique(mstrow(~inan) + hhmmssmat2serial(s.data.Time(~inan,:)));
counts      =  accumarray(subs, 1)-1;
% Aggregate by mst row
res = cached(:,{'Id','Date'});
res.Nconsolidated = uint32(accumarray(fix(un),  counts, [nmst,1]));
end

% Sampling
function res = sample(s,cached, opt)

% Sampling params
if nargin < 3
    opt.grid    = (9.5/24:5/(60*24):16/24)';
    opt.writeto = '.\data\TAQ\sampled\5min';
    opt.fmtname = 'S5m_%04d.mat';
end

nfile  = cached{end};
cached = cached{1};
% Number of observations per day
nobs   = double(s.mst.To - s.mst.From + 1);
nmst   = size(s.mst,1);

% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,iprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
iprice     = iprice ~= 1;
inan       = inan | iprice;

% STEP 3) Bad series/days
inan = inan | RunLength(cached.Isbadday,nobs);

% STEP 4) Filter out days with < 30min avg timestep or securities with 50% fewtrades days
inan = inan | RunLength(cached.Timestep,nobs);

% STEP 5) Clean prices - carried out in the median consolidation
% s.data.Price(inan) = NaN;

if ~all(inan)
    
    % STEP 6) Median prices for same timestamps
    mstrow           = RunLength((1:nmst)',nobs);
    [unTimes,~,subs] = unique(mstrow(~inan) + hhmmssmat2serial(s.data.Time(~inan,:)));
    price            = accumarray(subs, s.data.Price(~inan),[],@fast_median);
    
    % STEP 7) Sample on fixed grid (easier to match sp500)
    ngrid          = numel(opt.grid);
    [price, dates] = fixedsampling(unTimes, price, opt.grid(:));
    idx            = fix(dates);
    dates          = yyyymmdd2serial(double(s.mst.Date(idx))) + rem(dates,1);
    
    % Reorganize output
    res        = [];
    s.data     = table(dates,price,'VariableNames',{'Datetime','Price'});
    imst       = RunLength(idx);
    s.mst      = [s.mst(imst,'Id'), cached(imst, 'UnID'), s.mst(imst,'Date')];
    s.mst.From = (1:ngrid:size(s.data,1))';
    s.mst.To   = (ngrid:ngrid:size(s.data,1))';
    
    % Save
    fname = fullfile(opt.writeto,sprintf(opt.fmtname,nfile));
    save(fname, '-struct','s')
end
end

function res = betacomponents(s,cached)
spdays = cached{2};
spyret = cached{1};
useon  = numel(cached) == 4;
if useon
    onret = cached{3};
end
    
ngrid  = size(spyret{1},1);

% Dates and returns
dates  = s.data.Datetime;
ret    = [NaN; s.data.Price(2:end)./s.data.Price(1:end-1)-1];
idx    = [false; diff(rem(dates,1)) >= 0];
if useon
    [~,pos]   = ismember(s.mst(:,{'UnID','Date'}), onret(:,{'UnID','Date'}));
    ret(~idx) = onret.Onret(pos);
else
    % Keep all except overnight
    ret = ret(idx,:);
end
% Use a NaN when we don't have SPY returns
spyret = [NaN(ngrid,1); spyret];
days   = yyyymmdd2serial(double(s.mst.Date));
pos    = ismembc2(days, spdays) + 1;

% Map SP500 rets to stock rets
spret   = cat(1,spyret{pos});
prodret = spret.*ret;
subsID  = reshape(repmat(1:size(s.mst,1),ngrid,1),[],1);
ikeep   = ~isnan(prodret);

% Store results
res     = s.mst(:,{'Id','UnID','Date'});
res.Num = accumarray(subsID(ikeep), prodret(ikeep),[],[],NaN);
res.Den = accumarray(subsID(ikeep), spret(ikeep).^2,[],[],NaN);
end

function res = selrulecounts(s,cached)
nfile  = uint16(cached{end});
vnames = {'Date','Val','Count'};

% Sort mst
if ~issorted(s.mst.From)
    s.mst = sortrows(s.mst,'From');
end
dates = RunLength(s.mst.Date, double(s.mst.To-s.mst.From+1));

% G127
[g127,~,subs] = unique([dates,s.data.G127_Correction(:,1)],'rows');
g127          = table(g127(:,1),g127(:,2),accumarray(subs,1), 'VariableNames', vnames);

% Correction
[correction,~,subs] = unique([dates, s.data.G127_Correction(:,2)],'rows');
correction          = table(correction(:,1),correction(:,2),accumarray(subs,1), 'VariableNames', vnames);

% Condition
[condition,~,subs] = unique(table(dates, s.data.Condition));
condition.Count    = accumarray(subs,1);
condition.Properties.VariableNames = vnames;

% Null price
bins               = double(s.data.Price ~= 0);
[nullprice,~,subs] = unique([dates, bins],'rows');
nullprice          = table(nullprice(:,1), nullprice(:,2), accumarray(subs,1),'VariableNames', vnames);

% Null size
bins              = double(s.data.Volume ~= 0);
[nullsize,~,subs] = unique([dates, bins],'rows');
nullsize          = table(nullsize(:,1), nullsize(:,2), accumarray(subs,1),'VariableNames', vnames);

res = {g127, correction, condition, nullprice, nullsize, nfile};
end

function res = return_overnight(s,cached)
cached = cached{1};
% Calculate daily return with daily initial NaNs offset
pnan         = find(isnan(s.data.Price));
offset       = histc(pnan,[s.mst.From, s.mst.To+1]');
offset       = offset(1:2:end);
to           = s.mst.To;
from         = s.mst.From + uint32(offset);
cached.Ocret = s.data.Price(to)./s.data.Price(from)-1;
cached.Onret = (1 + cached.Totret)./(1 + cached.Ocret) - 1; 
res          = cached;
end

function res = countnullrets(s,cached)
res          = s.mst(:,{'Id','UnID','Date'});
idx          = mcolon(s.mst.From,1,s.mst.To);
subs         = RunLength(1:size(s.mst,1),s.mst.To-s.mst.From+1);
res.Nullrets = accumarray(subs(:), s.data.Price(idx),[],@(x) nnz(x(2:end)./x(1:end-1)==1));
end

%% Deprecated

% Check how many bad prices
function res = cleanprices(s,cached)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan   = selecttrades(s.data);
selper = nnz(inan)/numel(s.data.Price)*100;
% fprintf('Selection: %.1f\n',selper)

% STEP 2) Bad prices prices are 1.5x or .65x the daily median (net of selection nans)
[~,igoodprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
inan           = inan | igoodprice ~= 1;
% fprintf('Bad prices: %.1f\n',nnz(inan)/numel(s.data.Price)*100)

% STEP 3) Bad series/days
inan      = inan | RunLength(cached.Baddays,nobs);
badserper = nnz(inan)/numel(s.data.Price)*100;
% fprintf('Bad series/days: %.1f\n',badserper)

% STEP 4) Filter out days with < 30min avg timestep or securities with 50% fewtrades days
inan = inan | RunLength(cached.Timestep,nobs);
% fprintf('Timestep: %.1f\n',nnz(inan)/numel(s.data.Price)*100)

% STEP 5) Clean prices
s.data.Price(inan) = NaN;

res = NaN(1,2,'single');
if ~all(inan)
    % Returns
    ret  = s.data.Price(2:end)./s.data.Price(1:end-1)-1;
    % Excluding overnight
    ret(setdiff(s.mst.From-1,[-1,0])) = NaN;
    % MinMax
    [res(1),res(2)] = MinMaxElem(ret);
end
fprintf('\t\t\t\t\t\t\t %.1f - %.1f\n',selper,badserper)
end

function res = dailystats(s,cached)
cachedmst = cached{1};

% STEP 1) Selection
inan = selecttrades(s.data);

% Select non nans
prices = s.data.Price(~inan);
% Set to NaN overnight
[n, subs] = histc(find(~inan), [1; s.mst.To+1]);

% Daily median price
fill = NaN('single');
nmst = size(s.mst,1);
Med  = accumarray(subs, prices,[nmst,1], @fast_median, fill);

% Min max return
ret              = prices(2:end)./prices(1:end-1)-1;
nnovernight      = true(size(ret));
idx              = setdiff(cumsum(n(1:end-2)),0);
nnovernight(idx) = false;
ret              = ret(nnovernight);
subs             = subs(nnovernight);
% Collect results (note, nrets can be negative if it was only one price, which was selected out)
res          = cachedmst(:,{'Id','Date'});
res.Min      = accumarray(subs, ret,[nmst,1],@min,fill);
res.Max      = accumarray(subs, ret,[nmst,1],@max,fill);
res.MedPrice = Med;
% res.Nrets    = n(1:end-1)-1;
end

function res = avgtimestep(s,cached)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(s.mst,1);
% STEP 1) Selection
inan = selecttrades(s.data);

% STEP 2) Bad prices are < than .65x daily median or > than 1.5x daily median
[~,ibadprice] = histc(s.data.Price./RunLength(cached.MedPrice,nobs), [.65,1.51]);
ibadprice     = ibadprice ~= 1;
inan          = inan | ibadprice;

% STEP 3) Take median price if same timestamp
mstrow = RunLength((1:nmst)',nobs);
mstrow = mstrow(~inan);
times  = unique(mstrow + hhmmssmat2serial(s.data.Time(~inan,:)));

% Calculate average intraday unique time step after selection, bad prices and consolidation
subs         = uint32(fix(times));
l            = ones('double');
n            = accumarray(subs,      l, [nmst,1],   [], l);
openTime     = accumarray(subs,  times, [nmst,1], @min, l);
closeTime    = accumarray(subs,  times, [nmst,1], @max, l);
% openTime     = 9.5/24;
% closeTime    = 16/24;
res          = cached(:,{'Id','Date'});
res.Timestep = (closeTime - openTime)./(n-1);
end

function res = betas(s,cached)
spdays = cached{1}{2};
spyret = cached{1}{1};
ngrid  = size(spyret{1},1);

% Dates and returns
dates  = s.data.Datetime;
ret    = s.data.Price(2:end)./s.data.Price(1:end-1)-1;

% Keep all except overnight
idx    = diff(rem(dates,1)) >= 0;
ret    = ret(idx,:);

% Use a NaN when we don't have SPY returns
spyret = [NaN(ngrid,1); spyret];
days   = yyyymmdd2serial(double(s.mst.Date));
pos    = ismembc2(days, spdays) + 1;

% Map SP500 rets to stock rets
spret   = cat(1,spyret{pos});
prodret = spret.*ret;
subsID  = reshape(repmat(1:size(s.mst,1),ngrid,1),[],1);
ikeep   = ~isnan(prodret);
beta    = accumarray(subsID(ikeep), prodret(ikeep))./accumarray(subsID(ikeep), spret(ikeep).^2,[],[],NaN);

% Store results
res = s.mst(:,{'Id','UnID','Date'});
res.Beta = single(beta);
end