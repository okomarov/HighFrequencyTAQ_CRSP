function [stratret, tbstats, tbarets] = momstrat(tb)

% Unstack net returns
netRet = unstack(tb(:,{'UnID','Date','Netret'}),'Netret','UnID');
netRet = sortrows(netRet,'Date');
dates  = netRet.Date;
netRet = table2array(netRet(:,2:end));

% Unstack tot returns
totRet = unstack(tb(:,{'UnID','Date','Dayret'}),'Dayret','UnID');
totRet = sortrows(totRet,'Date');
totRet = table2array(totRet(:,2:end));

% Group by month
[~, ~, subs] = unique(dates./100);
perbalance = find([false; diff(subs)>0]);

N        = numel(dates);
stratret = NaN(N, 2);
for ii = 1:numel(perbalance)
    [netIshort, netIlong, totIshort,totIlong] = deal(false(1,size(totRet,2)));
    % Alive in the future
    ifut    = subs == ii+1;
    isalive = any(netRet(ifut,:));
    
    % Calculate past performance
    ipast    = subs == ii;
    
    % Net
    pastperf = netRet(ipast,isalive);
    pastperf(isnan(pastperf)) = 0;
    pastperf = prod(1 + pastperf)-1;
    netPtiles = prctile(pastperf,[10,90]);
    netIshort(isalive) = pastperf <= netPtiles(1);
    netIlong (isalive) = pastperf >= netPtiles(2);
    
    % Tot
    pastperf = totRet(ipast,isalive);
    pastperf(isnan(pastperf)) = 0;
    pastperf = prod(1 + pastperf)-1;
    totPtiles = prctile(pastperf,[10,90]);
    totIshort(isalive) = pastperf <= totPtiles(1);
    totIlong (isalive) = pastperf >= totPtiles(2);
    
    % Daily strat
    stratret(ifut,1) = nanmean([totRet(ifut,netIlong), -totRet(ifut,netIshort)],2);
    stratret(ifut,2) = nanmean([totRet(ifut,totIlong), -totRet(ifut,totIshort)],2);
end

plotdates = datetime(yyyymmdd2serial(dates),'ConvertFrom','datenum');
from      = find(~all(isnan(stratret),2),1,'first');
stratlvl  = [NaN(from-2,2); cumprod([ones(1,2); stratret(from:end,:)+1])];

if nargout == 0
    plot(plotdates, stratlvl)
    legend({'Net','Tot'})
end

[tbstats, tbarets] = stratstats(stratlvl(from-1:end,:),plotdates(from-1:end));

end

function [tbstats, tbarets] = stratstats(lvl,dates)
monthrets     = level2mrets(lvl,dates);
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