function [res, filename] = AnalyzeImom(fun, varnames, cached, path2data, debug, poolcores, varargin)
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
%   ANALYZE(..., PATH2DATA) If you wanna use other than '.\data\TAQ\'
%                           files (default), then specify a different
%                           PATH2DATA, e.g '.\data\TAQ\sampled\5min'
%
%   ANALYZE(..., CACHED) Some FUN might require pre-cached results which where
%                        run on the whole database.
%                        Check the specific sub-function for the format of the
%                        needed CACHED results.
%   ANALYZE(..., DEBUG) Run execution sequentially, i.e. not in parallel, to be
%                       able to step through the code in debug mode.
if nargin < 2 || isempty(varnames),  varnames  = '';             end
if nargin < 3,                       cached    = [];             end
if nargin < 4 || isempty(path2data); path2data = '.\data\TAQ\';  end
if nargin < 5 || isempty(debug);     debug     = false;          end
if nargin < 6 || isempty(poolcores); poolcores = 8;              end

fhandles = {@avgObsPerBucket
            @isEnoughObs
            @halfHourRet
            @halfHourVol};

[hasFunc, pos] = ismember(fun, cellfun(@func2str,fhandles,'un',0));
if ~hasFunc
    error('Unrecognized function "%s".', fun)
end
fun             = fhandles{pos};
projectpath     = fileparts(mfilename('fullpath'));
[res, filename] = blockprocess(fun,projectpath, varnames, cached,path2data,debug,poolcores,varargin{:});
end

%% Subfunctions
function res = avgObsPerBucket(s, cached, opt)
data  = s.res;

% Group by date - HHMMSS
key = uint64(data.Date)*1e6 + uint64(data.HHMMSS);

[~,pos, subs] = unique(key);
res = data(pos,{'Date','HHMMSS'});

res.Sum = accumarray(subs, uint64(data.Counts));
res.N   = accumarray(subs, data.Counts~=0);
end

function res = isEnoughObs(s, cached, opt)
data  = s.res;
ikeep = data.HHMMSS < opt.Time;
data  = data(ikeep,:);

% Group by id - date
key = int64(data.Id)*1e8 + int64(data.Date);

% Preallocate
idx = [true; logical(diff(key))];
res = data(idx,{'Id','Date'});

% Sum counts up to time
subs       = cumsum(idx);
res.Counts = accumarray(subs,data.Counts);

res.HasEnoughObs = res.Counts >= opt.Minobs;
end

function res = halfHourRet(s, cached, opt)
cached = cached{1};

% Filter permnos of interest
[idx,pos] = ismembIdDate(s.mst.Permno, s.mst.Date, cached.Permno, cached.Date);
s.mst     = s.mst(idx,:);
pdata     = mcolonint(s.mst.From, s.mst.To);
s.data    = s.data(pdata,:);
cached    = cached(pos(idx),:);

% Filter half-hour of interest
hhmmss = serial2hhmmss(s.data.Datetime);
idx    = in(hhmmss, opt.HalfHourRange);
s.data = s.data(idx,:);

% Reshape data
nobs    = 7;
nseries = size(s.mst,1);
price   = reshape(s.data.Price,nobs,nseries);

% Replace first or last price
if opt.HalfHourRange(1) == 93000
    ret = price(end,:)'./cached.FirstPrice-1;
elseif opt.HalfHourRange(2) == 160000
    ret = cached.LastPrice./price(1,:)'-1;
else
    ret = (price(end,:)./price(1,:))'-1;
end

% Store values
res         = s.mst(:,{'Date','Permno'});
fname       = sprintf('T%d', opt.HalfHourRange(1));
res.(fname) = ret;
end

function res = halfHourVol(s, cached, opt)
cached = cached{1};

% Filter permnos of interest
[idx,pos] = ismembIdDate(s.mst.Permno, s.mst.Date, cached.Permno, cached.Date);
s.mst     = s.mst(idx,:);
pdata     = mcolonint(s.mst.From, s.mst.To);
s.data    = s.data(pdata,:);
cached    = cached(pos(idx),:);

% Filter half-hour of interest
hhmmss = serial2hhmmss(s.data.Datetime);
idx    = in(hhmmss, opt.HalfHourRange);
s.data = s.data(idx,:);

% Reshape data
nobs    = 7;
nseries = size(s.mst,1);
price   = reshape(s.data.Price,nobs,nseries);

% Replace first or last price
if opt.HalfHourRange(1) == 93000
    row        = max(sum(isnan(price)),1);
    pos        = sub2ind(size(price),row, 1:nseries);
    price(pos) = cached.FirstPrice;
elseif opt.HalfHourRange(2) == 160000
    price(end,:) = cached.LastPrice;
end
ret = price(2:end,:)./price(1:end-1,:)-1;

res         = cached(:,{'Date','Permno'});
fname       = sprintf('RV5_%d', opt.HalfHourRange(1));

% Calculate vol
res.(fname) = nansum(ret.*ret)';
end