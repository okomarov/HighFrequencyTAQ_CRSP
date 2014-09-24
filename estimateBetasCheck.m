%% Check 1-day betas
% Calculate Betas with new interface
betas = estimateBetas(1,5,true);
betas = sortrows(betas,{'UnID','Date'});

% Check against previously calculated ones
load .\results\20140813_1818_betas.mat
res = sortrows(res,{'UnID','Date'});

% Visualize difference
abdiff = abs(single(betas.Beta)-res.Beta);
counts = histc(abdiff, eps*10.^(0:15));
bar((0:15)', counts)

%% Check if multiple days extends

% Estimate
lookback = 20;
betas = estimateBetas(lookback,30,false);

% Load spyders
spy = loadresults('spysampled30m');

% Pick one UnID at random
unids      = randsample(unique(betas.UnID),1);
path2data  = '.\data\TAQ\sampled\30min\';
master     = load(fullfile(path2data,'master'),'-mat','mst');
idx        = ismember(master.mst.UnID,unids);
master.mst = master.mst(idx,:);
data       = getTaqData(master.mst, [],[],[],[],path2data);

% Add spyders
[idx,pos] = ismember(data.Datetime, spy.Datetime);
data.Spy(idx,1) = spy.Price(pos(idx));

% Get rid of NaNs
data = data(~isnan(data.Price) & ~isnan(data.Spy),:);

% Index by day
data.Date   = fix(data.Datetime);
manual      = table(unique(data.Date),'VariableNames',{'Date'});
N           = size(manual,1);
manual.Beta = NaN(N,1);

% Convert to return
data.Ret    = NaN(size(data,1),1);
data.Spyret = NaN(size(data,1),1);
for ii = 1:N
    idx = data.Date == manual.Date(ii);
    prices = data.Price(idx);
    spy    = data.Spy(idx);
    data.Ret(idx)    = [NaN; prices(2:end)./prices(1:end-1)-1];
    data.Spyret(idx) = [NaN; spy(2:end)./spy(1:end-1)-1];
end

% Drop overnight NaNs
data = data(~isnan(data.Ret),:);

% Calculate Betas
for ii = lookback:N
    idx = ismembc(data.Date, manual.Date(ii-lookback+1:ii));
    num = data.Ret(idx)'*data.Spyret(idx);
    den = data.Spyret(idx)'*data.Spyret(idx);
    manual.Beta(ii) = num./den;
end
manual.Date = serial2yyyymmdd(manual.Date);

% Compare
compare = betas(betas.UnID == unids & ismember(betas.Date, manual.Date),:);




