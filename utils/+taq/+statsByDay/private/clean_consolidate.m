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
    medprice = RunLength(cached.MedPrice,nobs);
    iexclude = iexclude | ...
               s.data.Price >= medprice .* multiplier |...
               s.data.Price <= medprice ./ multiplier;
end

if opts.ExcludeBadDays
    iexclude = iexclude | RunLength(cached.Isbadday, nobs);
end

if all(iexclude)
    data = [];
    return
end

s.data.Id(:,1)     = RunLength(s.index.Id, nobs);
s.data.Permno(:,1) = RunLength(s.index.Permno, nobs);
s.data.Date(:,1)   = RunLength(s.index.Date, nobs);
s.data             = s.data(~iexclude,:);

if ~isempty(opts.ConsolidateTimestampType)
    data = consolidateTimestamp(s.data, opts.ConsolidateTimestampType);
end
end
