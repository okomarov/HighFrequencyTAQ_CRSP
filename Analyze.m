function [res, filename] = Analyze(fun, varnames, cached, path2data, debug, varargin)
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

fhandles = {@medianprice
        @badprices
        @ibadprices
        @consolidationcounts
        @NumTimeBuckets
        @sample
        @sampleFirstLast
        @VWAP
        @sampleSpy};
    
[hasFunc, pos] = ismember(fun, cellfun(@func2str,fhandles,'un',0));
if ~hasFunc
    error('Unrecognized function "%s".', fun)
end
fun             = fhandles{pos};
projectpath     = fileparts(mfilename('fullpath'));
[res, filename] = blockprocess(fun ,projectpath, varnames, cached,path2data,debug, varargin{:});
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
function res = badprices(s, cached, multiplier)
if nargin < 3, multiplier = []; end

cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Selection
inan = isInvalidTrade(s.data);

% STEP 2) Bad prices are x times bigger or smaller than daily median
% NOTE: this is a lookahead filter which in real time should be substituted
% with e.g. previous day close adjusted for company events
if ~isempty(multiplier)
    medprice   = RunLength(cached.MedPrice,nobs);
    igoodprice = s.data.Price./medprice < multiplier &...
                 medprice./s.data.Price < multiplier;
else
    igoodprice = true(size(inan));
end

% STEP 3) Bad days
res         = cached(:,{'Id','Date'});
subs        = uint32(RunLength((1:size(cached,1))',nobs));
res.Nbadsel = uint32(accumarray(subs,  inan));
res.Nbadtot = uint32(accumarray(subs,  inan | ~igoodprice));
end

function ibad = ibadprices(s, cached, multiplier)
% Service function to flag bad prices in Analyze
% 
% Inputs:
%   s.data - data table
%   s.mst  - master records table to data
%   s.idx  - cell array of tickers
%  
%   cached - table aligned to s.mst with 'MedPrice' and 'Isbadday' fields
%
%   edges  - double with [lb, ub]

if nargin < 3, multiplier = []; end

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1) Select irregular trades
ibad = isInvalidTrade(s.data);

% STEP 2) Select bad prices (far from daily median)
if ~isempty(multiplier)
    medprice   = RunLength(cached.MedPrice,nobs);
    igoodprice = s.data.Price./medprice < multiplier &...
                 medprice./s.data.Price < multiplier;
    ibad       = ibad | ~igoodprice;
end
end

% Count how many observations we loose in consolidation step
function res = consolidationcounts(s,cached,edges)
cached = cached{1};

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s,cached,edges);

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

function res = NumTimeBuckets(s,cached,edges)
cached = cached{1};
% Note: last bin is lb <= x <= ub since data ends at 16:00
grid   = ([9.5:0.5:15.5, 16.5])/24;
fun    = @(x) nnz(histcounts(x, grid));

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s, cached,edges);

% STEP 4) Number of time buckets that have a trade
nmst                  = size(s.mst,1);
mstrow                = RunLength((1:nmst)',nobs);
Time                  = hhmmssmat2serial(s.data.Time(~ibad,:));
cached.NumTimeBuckets = uint8(accumarray(mstrow(~ibad),Time,[nmst,1], fun));

res = cached(:,{'Id','Date','NumTimeBuckets'});
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
    s.mst.From = (1:ngrid:size(s.data,1))';
    s.mst.To   = (ngrid:ngrid:size(s.data,1))';
    
    % Save
    fname = fullfile(opt.writeto,sprintf(opt.fmtname,nfile));
    save(fname, '-struct','s')
end
end

% Sampling
function res = sampleFirstLast(s,cached,opt)
nfile  = cached{end};
cached = cached{1};
res    = [];

[price, times] = samplePrepare_(s,cached,opt);

if ~isempty(price)

    res     = cached(:,{'Date','Permno'});
    posdata = fix(times);
    
    % First
    idx                            = [true; logical(diff(fix(times)))];
    res.FirstPrice(posdata(idx),1) = price(idx); 
    res.FirstTime(posdata(idx),1)  = uint32(serial2hhmmss(times(idx))); 
    % Last
    idx                            = [idx(2:end); true];
    res.LastPrice(posdata(idx),1)  = price(idx); 
    res.LastTime(posdata(idx),1)   = uint32(serial2hhmmss(times(idx))); 
    
    res.File = repmat(uint16(nfile), size(res,1),1);
end
end
% Sampling
function res = VWAP(s,cached,opt)
% Volume weighted average price
nfile  = cached{end};
cached = cached{1};
res    = [];

[price, times, vol] = samplePrepare_(s,cached,opt);

if ~isempty(price)
    res = cached(:,{'Date','Permno'});
   
    % Accumulation subs 
    row    = fix(times);
    hhmmss = serial2hhmmss(times);
    col    = zeros(size(row));
    nedges = size(opt.edgesVWAP,1);
    for r = 1:nedges
        idx      = in(hhmmss, opt.edgesVWAP(r,:));
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

if isfield(opt,'edges')
    edges = opt.edges;
else
    edges = [];
end

% Number of observations per day
nobs = double(s.mst.To - s.mst.From + 1);

% STEP 1-3) Bad prices
ibad = ibadprices(s, cached, edges);

% STEP 4) Filter out days with < 30min avg timestep or securities with 50% fewtrades days
ibad = ibad | RunLength(cached.Isfewobs,nobs);

% STEP 5) Clean prices - carried out in the median consolidation
% s.data.Price(inan) = NaN;

% STEP 6) Volume-weighted average price for same timestamps
if ~all(ibad)
    nmst           = size(s.mst,1);
    mstrow         = RunLength((1:nmst)',nobs);
    [times,~,subs] = unique(mstrow(~ibad) + hhmmssmat2serial(s.data.Time(~ibad,:)));
    vol            = double(s.data.Volume)/100;
    voltot         = accumarray(subs, vol(~ibad));
    prices         = accumarray(subs, s.data.Price(~ibad).*vol(~ibad)) ./ voltot;
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

%% Deprecated

% % Identify bad prices
% function res = badprices(s, cached, edges)
% if nargin < 4, edges = []; end
% 
% cached = cached{1};
% 
% % Number of observations per day
% nobs = double(s.mst.To - s.mst.From + 1);
% 
% % STEP 1) Selection
% inan = isInvalidTrade(s.data);
% 
% % STEP 2) Bad prices are < than .5x daily median or > than 1.5x daily median
% if ~isempty(edges)
%     igoodprice = in(s.data.Price./RunLength(cached.MedPrice,nobs), edges);
% else
%     igoodprice = true(size(inan));
% end
% 
% % STEP 3) Bad days
% res         = cached(:,{'Id','Date'});
% subs        = uint32(RunLength((1:size(cached,1))',nobs));
% res.Nbadsel = uint32(accumarray(subs,  inan));
% res.Nbadtot = uint32(accumarray(subs,  inan | ~igoodprice));
% end