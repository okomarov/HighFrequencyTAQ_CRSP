function estimateHighLowSpread(LAG)
% beta = E[sum_j(log(high_j/low_j)^2 for j = 0:1)] and expectation as average over past LAG days
% Reference: 2012 Corwin, Schults - A Simple Way to Estimate Bid-Ask Spreads from Daily High and Low Prices - JF

% Load high low prices and cache
hl      = loadresults('estimateHighLowPrice','..\results');
hl      = hl(~isnan(hl.High) & ~isnan(hl.Low),:);
hl.High = double(hl.High);
hl.Low  = double(hl.Low);
gvar    = findgroups(hl.Permno);
rows    = accumarray(gvar,1);
hl      = cache2cell(hl, findgroups(hl.Permno));

for ii = 1:numel(hl)
    if rows(ii) < OPT_LAG_RUN+OPT_LAGDAY
        hl{ii}.Beta = NaN(rows(ii),1);
    else
        hl{ii}.Beta        = cumsum(log(hl{ii}.High ./ hl{ii}.Low).^2,1);
        hl{ii}.Beta(3:end) = hl{ii}.Beta(3:end) - hl{ii}.Beta(1:end-2);
        hl{ii}.Beta(1) = NaN;
        hl{ii}.Beta(2:end) = movmean(hl{ii}.Beta(2:end),[LAG-1,0]);
        hl{ii}.Beta(2:end) = hl{ii}.Beta(1:end-1);
        hl{ii}.Beta(1:LAG+1) = NaN;

        log(max(hl{ii}.High(1:end-1), hl{ii}.High(2:end)) ./ min(hl{ii}.Low(1:end-1), hl{ii}.Low(2:end))).^2
    end
end
end
