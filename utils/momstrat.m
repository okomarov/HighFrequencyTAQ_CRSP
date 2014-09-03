function stratret = momstrat(tb)

% Unstack net returns
netRet = unstack(tb(:,{'UnID','Date','Netret'}),'Netret','UnID');
netRet = sortrows(netRet,'Date');
dates  = netRet.Date;
netRet = table2array(netRet(:,2:end))+1;

% Unstack tot returns
totRet = unstack(tb(:,{'UnID','Date','Dayret'}),'Dayret','UnID');
totRet = sortrows(totRet,'Date');
totRet = table2array(totRet(:,2:end));

% Group by month
[~, ~, subs] = unique(dates./100);
prebalance = find([false; diff(subs)>0]);

N        = numel(dates);
stratret = NaN(N, 2);
for ii = 1:numel(prebalance)
    % Alive in the future
    ifut    = subs == ii+1;
    isalive = any(netRet(ifut,:));
    
    % Calculate past performance
    ipast    = subs == ii;
    
    % Net
    pastperf = netRet(ipast,isalive);
    pastperf(isnan(pastperf)) = 0;
    netPtiles = prctile(prod(1 + pastperf)-1,[10,90]);
    netIshort(isalive) = netRet(prebalance(ii),isalive) <= netPtiles(1);
    netIlong (isalive) = netRet(prebalance(ii),isalive) >= netPtiles(2); 
        
    % Tot
    pastperf = totRet(ipast,isalive);
    pastperf(isnan(pastperf)) = 0;
    totPtiles = prctile(prod(1 + pastperf)-1,[10,90]);
    totIshort(isalive) = totRet(prebalance(ii),isalive) <= totPtiles(1);
    totIlong (isalive) = totRet(prebalance(ii),isalive) >= totPtiles(2); 
    
    % Daily strat
    stratret(ifut,1) = nanmean([totRet(ifut,netIlong), -totRet(ifut,netIshort)],2);
    stratret(ifut,2) = nanmean([totRet(ifut,totIlong), -totRet(ifut,totIshort)],2);
end
plotdates = datetime(yyyymmdd2serial(dates),'ConvertFrom','datenum');
plot(plotdates, cumprod([ones(1,2); stratret(2:end,:)+1]))
legend({'Net','Tot'})
end