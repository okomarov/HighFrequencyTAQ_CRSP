function data = clean_consolidate(s, cached, opts)
% CLEAN_CONSOLIDATE Cleans invalid/bad trades and consolidates timestamp

if ~isempty(cached) && ~isequal(s.index.Date, cached.Date)
    error('Index and cache are not aligned.')
end

if nargin < 3 || isempty(opts)
    opts = struct('ExcludeBadDays',false, 'BadPriceMultiplier', [],'ConsolidateTimestampType',[]);
end

nobs     = double(s.index.To - s.index.From + 1);
iexclude = isInvalidTrade(s.data);

% Exclude prices > multiplier or prices < 1/multiplier
if ~isempty(opts.BadPriceMultiplier)
    medprice = RunLength(cached.MedianPrice,nobs);
    iexclude = iexclude | ...
               s.data.Price >= medprice .* opts.BadPriceMultiplier |...
               s.data.Price <= medprice ./ opts.BadPriceMultiplier;
end

if opts.ExcludeBadDays
    iexclude = iexclude | RunLength(cached.Isbadday, nobs);
end

if all(iexclude) || opts.TerminateEarly
    data = [];
    return
end

% Add id and date to data
data = s.data;
data.Id(:,1)     = RunLength(s.index.Id    , nobs);
data.Permno(:,1) = RunLength(s.index.Permno, nobs);
data.Date(:,1)   = RunLength(s.index.Date  , nobs);
data             = data(~iexclude,:);

if ~isempty(opts.ConsolidateTimestampType)
    data = consolidateTimestamp(data, opts.ConsolidateTimestampType);
end
end
