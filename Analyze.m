function [res, filename] = Analyze(fun, varnames, cached, path2data, debug, poolcores, varargin)
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
if nargin < 4 || isempty(path2data); path2data = '.\data\TAQ\';         end
if nargin < 5 || isempty(debug);     debug     = false;                 end
if nargin < 6 || isempty(poolcores); poolcores = 4;                     end

fhandles = {@medianprice
    @countBadPrices
    @ibadprices
    @consolidationcounts
    @countPerTimeBucket
    @sample
    @sampleFirstLast
    @betacomponents_refresh
    @betacomponents_preav
    @VWAP
    @sampleSpy
    @testLoadSpeed
    @maxRecordsPerSec
    @countNonNullRet};

[hasFunc, pos] = ismember(fun, cellfun(@func2str,fhandles,'un',0));
if ~hasFunc
    error('Unrecognized function "%s".', fun)
end
fun             = fhandles{pos};
projectpath     = fileparts(mfilename('fullpath'));
[res, filename] = blockprocess(fun ,projectpath, varnames, cached,path2data,debug,poolcores, varargin{:});
end
%% Subfunctions

% Median price (net of the first selection step)
function res = medianprice(s,cached)
cached = cached{1};

% STEP 1) Selection
inan = isInvalidTrade(s.data);

% Prepare for daily median
nobs = double(s.mst.To - s.mst.From + 1);
nmst = size(cached,1);
subs = uint32(RunLength((1:nmst)',nobs));

% Daily median price
cached.MedPrice      = accumarray(subs(~inan), s.data.Price(~inan),[nmst,1], @fast_median);
idx                  = cached.MedPrice == 0;
cached.MedPrice(idx) = NaN;

res = cached;
end

% Identify bad prices
function res = countBadPrices(s, cached, multiplier)
if nargin < 3, multiplier = []; end

cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

[ibad, invalid] = ibadprices(s, cached, multiplier);

% Counts
res         = cached(:,{'Id','Date'});
subs        = uint32(RunLength((1:size(cached,1))',nobs));
res.Ntot    = nobs;
res.Nbadsel = uint32(accumarray(subs,  invalid));
res.Nbadtot = uint32(accumarray(subs,  ibad));
end

function [ibad, invalid] = ibadprices(s, cached, multiplier)
% Service function to flag bad prices in Analyze
%
% Inputs:
%   s.data - data table
%   s.mst  - master records table to data
%   s.idx  - cell array of tickers
%
%   cached - table aligned to s.mst with 'MedPrice'
%
%   multiplier - prices are bad if x times the median

if nargin < 3, multiplier = []; end

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Select irregular trades
invalid = isInvalidTrade(s.data);

% STEP 2) Select bad prices (far from daily median)
if ~isempty(multiplier)
    medprice = RunLength(cached.MedPrice,nobs);
    ibad     = s.data.Price./medprice >= multiplier |...
               medprice./s.data.Price >= multiplier;
    ibad     = ibad | invalid;
else
    ibad = invalid;
end
end

% Count how many observations we loose in consolidation step
function res = consolidationcounts(s,cached,multiplier)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s,cached,multiplier);

% STEP 4) Count how many observations we loose from median consolidation
nmst              = size(s.mst,1);
mstrow            = RunLength((1:nmst)',nobs);
% Count by mst row and timestamp
[un,~,subs]       = unique(mstrow(~ibad) + hhmmssmat2serial(s.data.Time(~ibad,:)));
counts            = accumarray(subs, 1)-1;
% Aggregate by mst row
res               = cached(:,{'Id','Date'});
res.Nconsolidated = uint32(accumarray(fix(un),  counts, [nmst,1]));
end

function res = countPerTimeBucket(s,cached,opt)
% Counts number of trades within a time bucket
nfile  = cached{end};
cached = cached{1};

nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s, cached, opt.BadPriceMultiplier);

% Note: last bin is lb <= x <= ub since data ends at 16:00
nmst           = (1:size(s.mst,1))';
[grid, mstrow] = ndgrid(opt.Grid,nmst);
times          = RunLength(nmst, nobs) + hhmmssmat2serial(s.data.Time);
times          = times(~ibad);
times          = times([true; logical(diff(times))]);
counts         = histcounts(times, mstrow(:)+grid(:));

% Exclude day-to-day counts
n               = size(grid,1);
counts(n:n:end) = [];
mstrow(n:n:end) = [];
grid(n:n:end)   = [];

res        = cached(mstrow,{'Id','Date'});
res.HHMMSS = uint32(serial2hhmmss(grid));
res.Counts = uint32(counts(:));

fname = fullfile('.\data\TAQ\count',sprintf('ntrades30m_f%04d.mat',nfile));
save(fname, 'res')
res   = [];
end

% Sampling
function res = sample(s,cached,opt)
nfile  = cached{end};
cached = cached{1};

[price, times] = samplePrepare_(s,cached,opt);

if ~isempty(price)
    % STEP 7) Sample on fixed grid (easier to match sp500)
    ngrid          = numel(opt.grid);
    [price, dates] = fixedsampling(times, price, opt.grid(:));
    idx            = fix(dates);
    dates          = yyyymmdd2serial(double(s.mst.Date(idx))) + rem(dates,1);

    % Reorganize output
    res        = [];
    s.data     = table(dates,price,'VariableNames',{'Datetime','Price'});
    imst       = RunLength(idx);
    s.mst      = [s.mst(imst,'Id'), cached(imst, 'Permno'), s.mst(imst,'Date')];
    s.mst.From = uint32((1:ngrid:size(s.data,1))');
    s.mst.To   = uint32((ngrid:ngrid:size(s.data,1))');

    % Save
    fname = fullfile(opt.writeto,sprintf(opt.fmtname,nfile));
    save(fname, '-struct','s')
end
end

% Sampling
function res = sampleFirstLast(s,cached,multiplier)
nfile  = cached{end};
cached = cached{1};
res    = [];

nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s, cached, multiplier);

if ~all(ibad)
    res = cached(:,{'Date','Id'});

    % Indexing into mst
    nmst   = size(s.mst,1);
    mstrow = RunLength((1:nmst)',nobs);
    mstrow = mstrow(~ibad);
    time   = uint32(double(s.data.Time(~ibad,:))*[10000;100;1]);
    price  = s.data.Price(~ibad);

    % First
    idx                   = [true; logical(diff(mstrow))];
    pos                   = mstrow(idx);
    res.FirstPrice(pos,1) = price(idx);
    res.FirstTime(pos,1)  = time(idx);
    % Last
    idx                   = [idx(2:end); true];
    pos                   = mstrow(idx);
    res.LastPrice(pos,1)  = price(idx);
    res.LastTime(pos,1)   = time(idx);
end
end

function res = betacomponents_refresh(s,cached,opt)
spydates = cached{2};
spymst   = cached{3};
cached   = cached{1};
nmst     = size(s.mst,1);

% Get all prices
[price, times] = samplePrepare_(s,cached,opt);
price          = cache2cell([rem(times,1) double(price)],fix(times),[],[nmst,1]);
iempty         = cellfun('isempty',price);

spy = getSpy_(spymst,opt.TimestampConsolidation);

res = cached(:,{'Permno','Date'});
[Num, Den] = deal(NaN(nmst,1));

for ii = 1:nmst
    idx = s.mst.Date(ii) == spydates;
    if ~iempty(ii) && res.Permno(ii) ~= 0 && any(idx)

        [ia,ispy] = sample_refresh(price{ii}(:,1), spy{idx}(:,1),opt.MinRefreshStep);

        % returns
        ret    = diff(log(price{ii}(ia,2)));
        retspy = diff(log(spy{idx}(ispy,2)));

        if opt.HasOvernight
            ret    = [cached.RetCO(ii); ret];
            retspy = [spymst.RetCO(idx); retspy];
        end

        Num(ii) = sum(ret.*retspy);
        Den(ii) = sum(retspy.*retspy);
    end
end
res.Num = Num;
res.Den = Den;
end

function res = betacomponents_preav(s,cached,opt)
spydates = cached{2};
spymst   = cached{3};
cached   = cached{1};
nmst     = size(s.mst,1);

% Get all prices
[price, times] = samplePrepare_(s,cached,opt);
price          = cache2cell([rem(times,1) double(price)],fix(times),[],[nmst,1]);
iempty         = cellfun('isempty',price);

spy = getSpy_(spymst,opt.TimestampConsolidation);

res = cached(:,{'Permno','Date'});
[Num, Den] = deal(NaN(nmst,1));

for ii = 1:nmst
    idx = s.mst.Date(ii) == spydates;
    if ~iempty(ii) && res.Permno(ii) ~= 0 && any(idx)

        [pA,pS] = sample_preaverage(price{ii}(:,1), spy{idx}(:,1),...
                                    price{ii}(:,2), spy{idx}(:,2), opt.MinRefreshStep);

        % returns
        ret    = diff(log(pA));
        retspy = diff(log(pS));

        if opt.HasOvernight
            ret    = [cached.RetCO(ii); ret];
            retspy = [spymst.RetCO(idx); retspy];
        end

        Num(ii) = sum(ret.*retspy);
        Den(ii) = sum(retspy.*retspy);
    end
end
res.Num = Num;
res.Den = Den;
end

function spy = getSpy_(mst,consolidationType)
% Get spy prices

nspy = size(mst,1);
spy  = cell(nspy,1);
for ii = 1:nspy
    spyprice = getTaqData([],[],[],[],[],'..\data\TAQ\',mst(ii,:),false);

    % STEP 1) Select regular trades
    spyprice = spyprice(~isInvalidTrade(spyprice),:);

    % STEP 2) Median prices for same timestamps
    [spytime,~,subs] = unique(rem(spyprice.Datetime,1));

    switch consolidationType
        case 'volumeWeighted'
            vol      = double(spyprice.Volume)/100;
            voltot   = accumarray(subs, vol);
            spyprice = accumarray(subs, double(spyprice.Price).*vol) ./ voltot;
        case 'median'
            spyprice = accumarray(subs, double(spyprice.Price),[],@fast_median);
    end

    spy{ii} = [spytime,spyprice];
end
end

% Sampling
function res = VWAP(s,cached,opt)
% Volume weighted average price
nfile  = cached{end};
cached = cached{1};
res    = [];

if ~isfield(opt,'inclusion')
    opt.inclusion = '[)';
end

[price, times, vol] = samplePrepare_(s,cached,opt);

if ~isempty(price)
    res = cached(:,{'Date','Permno'});

    % Accumulation subs
    row    = fix(times);
    hhmmss = serial2hhmmss(times);
    col    = zeros(size(row));
    nedges = size(opt.edgesVWAP,1);
    for r = 1:nedges
        idx      = in(hhmmss, opt.edgesVWAP(r,:),opt.inclusion);
        col(idx) = r;
    end
    subs = [row, col];

    sz = [size(res,1), nedges];

    % Drop uncategorized
    ikeep = col ~= 0;
    VWAP  = accumarray(subs(ikeep,:), price(ikeep).*vol(ikeep), sz) ./...
            accumarray(subs(ikeep,:), vol(ikeep), sz);

    % Keep only record for which we have prices
    row       = unique(row);
    res       = res(row,:);
    VWAP      = VWAP(row,:);
    VWAPnames = arrayfun(@(x) sprintf('T%d',x),opt.edgesVWAP(:,1),'un',0);
    res       = [res, array2table(VWAP,'VariableNames',VWAPnames)];
    res.File  = repmat(uint16(nfile), size(res,1),1);
end
end

function [prices,times, voltot] = samplePrepare_(s,cached,opt)
if isfield(opt,'BadPriceMultiplier')
    multiplier = opt.BadPriceMultiplier;
else
    multiplier = [];
end

if isfield(opt,'TimestampConsolidation')
    consolidationType = opt.TimestampConsolidation;
else
    consolidationType = 'volumeWeighted';
end

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s, cached, multiplier);

% Bad days
ibad = ibad | RunLength(cached.Isbadday, nobs);

if ~all(ibad)
    % Timestamp consolidation of price
    nmst           = size(s.mst,1);
    mstrow         = RunLength((1:nmst)',nobs);
    [times,~,subs] = unique(mstrow(~ibad) + hhmmssmat2serial(s.data.Time(~ibad,:)));
    prices         = s.data.Price(~ibad);
    vol            = double(s.data.Volume(~ibad))/100;
    voltot         = accumarray(subs, vol);
    switch consolidationType
        case 'first'
            idx    = [true; logical(diff(subs))];
            prices = prices(idx);
        case 'median'
            prices = accumarray(subs, prices,[],@fast_median);
        case 'volumeWeighted'
            prices = accumarray(subs, prices.*vol) ./ voltot;
        case 'skip'
            % Do not consolidate price
    end
else
    prices = [];
    times  = [];
    voltot = [];
end
end

function res = sampleSpy(s,~,opt)

% Keep spy only
id = find(strcmpi(s.ids,'SPY'),1);
if isempty(id)
    [Datetime,Price] = deal(zeros(0,1));
    res              = table(Datetime,Price);
    return
end
idx   = s.mst.Id == id;
s.mst = s.mst(idx,:);
s.mst = sortrows(s.mst, {'Date','From'});

% Loop for each day
nmst = size(s.mst,1);
res  = cell(nmst,1);
for r = 1:size(s.mst,1);
    from = s.mst.From(r);
    to   = s.mst.To(r);
    data = s.data(from:to,:);

    % STEP 1) Select regular trades
    data = data(~isInvalidTrade(data),:);

    % STEP 2) Median prices for same timestamps
    [unTimes,~,subs] = unique(hhmmssmat2serial(data.Time));
    Price            = accumarray(subs, double(data.Price),[],@fast_median);

    % STEP 3) Sample on fixed grid
    [Price, dates] = fixedsampling(unTimes, Price, opt.grid(:));
    Datetime       = yyyymmdd2serial(s.mst.Date(r)) + dates;

    res{r} = table(Datetime,Price);
end
res = cat(1,res{:});
end

function res = maxRecordsPerSec(s,~,~)
nobs  = double(s.mst.To - s.mst.From + 1);
dates = RunLength(yyyymmdd2serial(s.mst.Date),nobs);

[unDatetimes,~,gDatetime] = unique(dates + hhmmssmat2serial(s.data.Time));
[unDates,~,gDate]         = unique(fix(unDatetimes));
Num                       = accumarray(gDate,accumarray(gDatetime, 1),[],@max);
res                       = table(unDates, Num,'VariableNames',{'Date','Records'});
end

function res = countNonNullRet(s,cached)
cached = cached{1};
if ~isequal(s.mst.From, cached.From)
    error('Unequal mst and cached.')
end
res = cached(:,{'Id','Permno','Date'});

ret           = [diff(s.data.Price); NaN];
ret(s.mst.To) = NaN;
nonzero       = ret ~= 0 & ~(isnan(ret));

nobs           = double(s.mst.To - s.mst.From + 1);
nmst           = size(s.mst,1);
mstrow         = RunLength((1:nmst)',nobs);
ikeep          = nonzero & ~isInvalidTrade(s.data);
res.NonNullRet = accumarray(mstrow(ikeep), 1,[nmst,1]);
end
%% Utility functions
function out = testLoadSpeed(varargin)
out = [];
end