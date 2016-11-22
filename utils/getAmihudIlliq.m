function [out,unId,unDt] = getAmihudIlliq(id, dt, prc, ret, vol, n)
% GETAMIHUDILLIQ Estimate the end-of-month (log) Amihud illiquidity measure
%
%   Ranks from liquid (smaller values) to illiquid (bigger values)
%   
%   N specifies the number of months of data to use for the lookback

prc = double(prc);
ret = double(ret);
vol = double(vol);

% Amihud illiq |r_t|/(price_t.*volume_t)
illiq               = abs(ret ./ (prc.* vol));
illiq(isinf(illiq)) = NaN;
[illiq, unId, unDt] = my_unstack(id,dt,illiq);

% End-of-month rolling average of a year of daily observations
[mdate, pos, midx] = unique(unDt/100,'last');
nmonths            = numel(mdate);

% Fast rolling average, equivalent to:
%
% out = NaN(nmonths, size(illiq,2));
% 
% for ii = n:nmonths
%     idx       = ismember(midx, ii-n+1:ii);
%     out(ii,:) = nanmean(illiq(idx,:),1);
% end
% out = log(out);

% a) monthly cumulative sums
[nobs, csum] = deal(zeros(nmonths, size(illiq,2)));
nonnan       = cache2cell(~isnan(illiq),midx);
illiq        = cache2cell(illiq,midx);
for ii = 1:nmonths
    csum(ii+1,:) = sum(illiq{ii},1)  + csum(ii,:);
    nobs(ii+1,:) = sum(nonnan{ii},1) + nobs(ii,:);
end
csum = csum(2:end,:);
nobs = nobs(2:end,:);

% b) Cut cumulative sum at n
csum(n+1:end,:) = csum(n+1:end,:) - csum(1:end-n,:);
nobs(n+1:end,:) = nobs(n+1:end,:) - nobs(1:end-n,:);

% c) Mean
out = csum./nobs;

out             = log(out);
out(isinf(out)) = NaN;
end

function [panel, unId, unDt] = my_unstack(id,dt,val)
[unId, trash, col] = unique(id);
[unDt, trash, row] = unique(dt);

m     = numel(unDt);
n     = numel(unId);
panel = zeros(m,n);

pos        = sub2ind([m,n],row,col);
panel(pos) = val(:);
end