function sampleDataCheck(idtype)
narginchk(1,1)
% Sample random permno and date
idtype = upperfirst(idtype);
master = load(fullfile('.\data\TAQ\sampled\5min\','master'),'-mat');
id     = randsample(master.mst.(idtype),1);
mst    = master.mst(ismember(master.mst.(idtype), id),:);
date   = randsample(mst.Date,1);

% Get unsampled
unsampled         = getTaqData(idtype,id,date,date);
unsampled         = unsampled(~selecttrades(unsampled),:);
[Datetime,~,subs] = unique(unsampled.Datetime);
Price             = accumarray(subs, unsampled.Price,[],@fast_median);

% Sample with MFE
pt                         = 'D:\TAQ\HFbetas\utils\MFE'; addpath(pt) 
freq                       = 5;
grid                       = (9.5/24:freq/(60*24):16/24)';
[MfePrice, ~, MfeDatetime] = realized_price_filter(double(Price), Datetime,'unit','fixed', grid+yyyymmdd2serial(date));

% Same sampled day
sampled = getTaqData(idtype,id,date,date,[],'.\data\TAQ\sampled\5min');

% Compare
tb                           = [sampled, table(MfeDatetime,MfePrice)];
tb.MfePrice(isnan(tb.Price)) = NaN;

plot(serial2datetime(tb.Datetime),tb{:,{'Price','MfePrice'}})
end