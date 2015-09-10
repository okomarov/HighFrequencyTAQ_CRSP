function [tbstats, tbarets] = stratstats(dates, ret, freq)
if nargin < 3 || isempty(freq)
    freq = 'd';
end
if ~isdatetime(dates)
    dates = yyyymmdd2datetime(dates);
end
switch freq
    case 'd'
        scale = 252;
    case 'm'
        scale = 12;
end
tbstats       = table();

sz  = size(ret);
lvl = ret2lvl(ret);

fun          = @(x) hac(ones(sz(1),1), ret(:,x),'intercept',false,'display','off');
[~,se,coeff] = arrayfun(fun,1:sz(2));

tbstats.Avgret    = coeff(:);
tbstats.Se        = se(:);
tbstats.Pval(:,1) = tcdf(-abs(coeff./se),sz(1)-1)*2; 
tbstats.Annret    = lvl(end,:)'.^(1/years(dates(end)-dates(1)))-1;
tbstats.Annstd    = nanstd(ret)'*sqrt(scale);
tbstats.Downstd   = nanstd(ret > 0 .* ret)' * sqrt(scale);
tbstats.Minret    = nanmin(ret)';
tbstats.Maxret    = nanmax(ret)';
tbstats.IR        = tbstats.Annret./tbstats.Annstd;
[mdd,imdd]        = maxdrawdown(lvl);
tbstats.Mdd       = mdd(:);
tbstats.Reclen    = days(dates(imdd(end,:))-dates(imdd(1,:)));
tbstats.Sortino   = tbstats.Annret./tbstats.Downstd;

if nargout == 2
    tbarets = level2arets(lvl,dates);
end
end

function lvl = ret2lvl(ret)
inan      = isnan(ret);
ret(inan) = 0;
from      = find(~inan,1,'first')-1;
c         = size(ret,2);
lvl       = [NaN(from-1,c); 
             cumprod([ones(1,c); ret(from+1:end,:)+1])];
end

function arets = level2arets(lvl,dates)
sz    = size(lvl);
rsub  = repmat(cumsum([1; logical(diff(year(dates)))]), 1, sz(2));
csub  = repmat(1:sz(2),sz(1),1);
arets = accumarray([rsub(:),csub(:)],lvl(:), [],@(x) x(end)./x(1)-1); 
end