function [stratret, tbstats, tbarets] = zeroptf(tb, forceplot)
% [stratret, tbstats, tbarets] = zeroptf(tb)
% 
%   TB should have 'UnID', 'Date', 'Score' and 'Ret' variables

if nargin < 2 || isempty(forceplot), forceplot = false; end

% Remove rows with no score
tb = tb(~isnan(tb.Score),:);

warning off MATLAB:table:ModifiedVarnames

% Unstack scores
score = unstack(tb(:,{'UnID','Date','Score'}),'Score','UnID');
score = sortrows(score,'Date');
dates = uint32(score.Date);
score = table2array(score(:,2:end));

% Unstack returns
ret = unstack(tb(:,{'UnID','Date','Ret'}),'Ret','UnID');
ret = sortrows(ret,'Date');
ret = table2array(ret(:,2:end));

warning on MATLAB:table:ModifiedVarnames

% Beginning of month rebalancing scheme
[~, ~, subs] = unique(dates./100);
rebdate      = find([false; diff(subs)>0]);

N        = numel(dates);
stratret = NaN(N, 1);
for ii = 1:numel(rebdate)
    [iShort, iLong] = deal(false(1,size(ret,2)));
    % Alive on rebalancing date
    ifut    = subs == ii+1;
    isalive = ~isnan(score(rebdate(ii),:));
    
    % Score ranking
    scores          = score(rebdate(ii)-1,isalive);
    ptiles          = prctile(scores,[10,90]);
    iShort(isalive) = scores <= ptiles(1);
    iLong (isalive) = scores >= ptiles(2);
       
    % Equal weighted strategy
    stratret(ifut) = nanmean([ret(ifut,iLong), -ret(ifut,iShort)],2);
end

plotdates = yyyymmdd2datetime(dates);
from      = find(~isnan(stratret),1,'first');
stratlvl  = [NaN(from-2,1); cumprod([1; stratret(from:end)+1])];

if nargout == 0 || forceplot
    plot(plotdates, stratlvl)
end

[tbstats, tbarets] = stratstats(stratlvl(from-1:end,:),plotdates(from-1:end));
end

function [tbstats, tbarets] = stratstats(lvl,dates)
monthrets          = level2mrets(lvl,dates);
n                  = numel(monthrets);
[~,se,coeff]       = hac(ones(n,1), monthrets,'intercept',false,'display','off');
tbstats.Monret     = coeff;
tbstats.Pval       = tcdf(-abs(coeff/se),n-1)*2; 
tbstats.Annret     = lvl(end,:)'.^(1/years(dates(end)-dates(1)))-1;
tbstats.Annstd     = std(monthrets)'*sqrt(12);
tbstats.Downstd    = std(monthrets > 0 .* monthrets)' * sqrt(12);
tbstats.Minret     = min(monthrets)';
tbstats.Maxret     = max(monthrets)';
tbstats.IR         = tbstats.Annret./tbstats.Annstd;
[tbstats.Mdd,imdd] = maxdrawdown(lvl);
tbstats.Mdd        = tbstats.Mdd';
tbstats.Reclen     = days(dates(imdd(end,:))-dates(imdd(1,:)));
tbstats.Sortino    = tbstats.Annret./tbstats.Downstd;

tbstats = struct2table(tbstats);
tbarets = table(unique(year(dates)), level2arets(lvl,dates),'VariableNames',{'Year','Ret'});
end

function mrets = level2mrets(lvl,dates)
sz    = size(lvl);
rsub  = repmat(cumsum([1; logical(diff(month(dates)))]), 1, sz(2));
csub  = repmat(1:sz(2),sz(1),1);
mrets = accumarray([rsub(:),csub(:)],lvl(:), [],@(x) x(end)./x(1)-1); 
end

function arets = level2arets(lvl,dates)
sz    = size(lvl);
rsub  = repmat(cumsum([1; logical(diff(year(dates)))]), 1, sz(2));
csub  = repmat(1:sz(2),sz(1),1);
arets = accumarray([rsub(:),csub(:)],lvl(:), [],@(x) x(end)./x(1)-1); 
end