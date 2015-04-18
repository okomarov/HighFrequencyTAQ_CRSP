function getSpyCheck(freq)

% Take last price from spy sampled
spy = getSpy(freq);
spy = sortrows(spy,'Datetime');
pos = [find(diff(int64(serial2yyyymmdd(spy.Datetime))) ~= 0); size(spy,1)];
spy = spy(pos,:);

% Load dsfquery daily prices
dsfquery = loadresults('dsfquery');
dsfquery = dsfquery(dsfquery.Permno == 84398,:);
dsfquery  = sortrows(dsfquery, 'Date');

% Intersect dates
[~,ia,ib] = intersect(serial2yyyymmdd(spy.Datetime), dsfquery.Date);
spy = spy(ia,:);
dsfquery = dsfquery(ib,:);

% Visual inspection
subplot(211)
plot(serial2datetime(spy.Datetime), spy.Price,...
     yyyymmdd2datetime(dsfquery.Date), abs(dsfquery.Prc))
legend({'Sampled','Dsfquery'})

p2r = @(x) x(2:end)./x(1:end-1)-1;
perATE = abs(p2r(spy.Price) - p2r(abs(dsfquery.Prc)));
subplot(212);
boxplot(perATE)
str = {'Absolute Tracking Error'
       sprintf('%-10s%5.4f','mean:',mean(perATE))
       sprintf('%-10s%5.4f','std:',std(perATE))};
title(str)
end
