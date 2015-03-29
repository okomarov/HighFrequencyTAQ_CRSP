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
if ~exist('master','var')
    master = load(fullfile(path2data,'master'),'-mat','mst');
end
idx  = ismember(master.mst.UnID,unids);
data = getTaqData(master.mst(idx,:), [],[],[],[],path2data);

% Add spyders
data.Spy        = NaN(size(data,1),1);
[idx,pos]       = ismember(data.Datetime, spy.Datetime);
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
    idx              = data.Date == manual.Date(ii);
    prices           = data.Price(idx);
    spyprices        = data.Spy(idx);
    data.Ret(idx)    = [NaN; prices(2:end)./prices(1:end-1)-1];
    data.Spyret(idx) = [NaN; spyprices(2:end)./spyprices(1:end-1)-1];
end

% Drop overnight NaNs and empty days
data   = data(~isnan(data.Ret),:);
manual = manual(ismembc(manual.Date, data.Date),:);

% Calculate Betas
N = size(manual,1);
for ii = lookback:N
    idx = ismembc(data.Date, manual.Date(ii-lookback+1:ii));
    num = data.Ret(idx)'*data.Spyret(idx);
    den = data.Spyret(idx)'*data.Spyret(idx);
    manual.Beta(ii) = num./den;
end
manual.Date = serial2yyyymmdd(manual.Date);

% Compare
compare = betas(betas.UnID == unids & ismember(betas.Date, manual.Date),:);

% Visual inspection
plot(yyyymmdd2datetime(compare.Date), compare.Beta,...
     yyyymmdd2datetime(manual.Date), manual.Beta)
legend({'Automated','Manual'})

abdiff = abs(compare.Beta - manual.Beta);
counts = histc(abdiff, eps*10.^(0:15));
bar((0:15)', counts)
set(gca,'XLim',[-1 16],'Xtick',0:15)