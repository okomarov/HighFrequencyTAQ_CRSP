function mst = selectAndFilterTrades()

% Load big master file
path2data = '.\data\TAQ';
load(fullfile(path2data,'master'),'-mat')

% Map unique ID to mst
testname = 'masterPermno';
try
    res = loadresults(testname);
catch
    res = mapPermno2master;
end
[~,pos]    = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Permno = res.Permno(pos);

% Median price
testname = 'medianprice';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date'}));
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.MedPrice = res.MedPrice(pos);

% Nobs
mst.Nobs = mst.To - mst.From +1;

% Bad prices days
testname = 'badprices';
try
    res = loadresults(testname);
catch
    dailycut = 0.5;
    res      = Analyze(testname,[],mst(:, {'File','Id','Date','MedPrice'}),[],[],dailycut);
end
[~,pos]      = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Nbadsel  = res.Nbadsel(pos);
mst.Nbadtot  = res.Nbadtot(pos);
mst.Isbadday = res.Isbadday(pos);

% Bad series
subs         = mst.Permno + 1; % Shift by one to avoid 0 permno problems
sz           = [max(mst.Permno)+1,1];
totbad       = accumarray(subs(mst.Isbadday), mst.Nobs(mst.Isbadday),sz);
totobs       = accumarray(subs, mst.Nobs, sz);
threshold    = 0.1;
badseries    = totbad./totobs > threshold;
mst.Isbadday = mst.Isbadday | badseries(subs);

% Count losing obs with timestamp consolidation
testname = 'consolidationcounts';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],mst(:, {'File','Id','Date','MedPrice','Isbadday'}));
end
[~,pos]           = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.Nconsolidated = res.Nconsolidated(pos);

% Count number of time buckets in a day that have a trade
testname = 'NumTimeBuckets';
try
    res = loadresults(testname);
catch
    res = Analyze('NumTimeBuckets',[],mst(:, {'File','Id','Date','MedPrice','Isbadday'}));
end
[~,pos]            = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.NumTimeBuckets = res.NumTimeBuckets(pos);

% Select with minimum number of observations
subs            = mst.Permno + 1;
mst.Ngoodtrades = mst.Nobs - mst.Nbadtot - mst.Nconsolidated;
mst.Isfewobsday = mst.NumTimeBuckets < 7 | mst.Ngoodtrades < 40;
fewobsdays      = accumarray(subs, mst.Isfewobsday);
totdays         = accumarray(subs, 1);
threshold       = 0.5;
badseries       = fewobsdays./totdays > threshold;
mst.Isfewobs    = mst.Isfewobsday | badseries(subs);

end