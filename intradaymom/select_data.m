OPT_EDGES_BAD_PRICES = [.5,1.5];

%% Import data

% Index data
datapath = '..\data\TAQ\sampled\5min\new';
master = load(fullfile(datapath,'master'),'-mat');

% Sample first and last price
try
    price_fl = loadresults('sampleFirstLast','..\results');
catch
    mst = selectAndFilterTrades(OPT_EDGES_BAD_PRICES);
    if isempty(OPT_EDGES_BAD_PRICES)
        mst = mst(:, {'File','Id','Permno','Date','Isbadday','Isfewobs'});
    else
        mst = mst(:, {'File','Id','Permno','Date','MedPrice', 'Isbadday','Isfewobs'});
    end
    price_fl = Analyze('sampleFirstLast',[],mst,[],[], struct('edges',[]));
end
[idx,pos] = ismembIdDate(price_fl.Permno, price_fl.Date, master.mst.Permno, master.mst.Date);
tmp(pos(idx),:) = price_fl(idx,:);
price_fl = tmp;
price_fl.File = master.mst.File;

% Common shares
shrcd      = loadresults('shrcd','../results');
idx        = iscommonshare(master.mst, shrcd);
master.mst = master.mst(idx,:);
idx        = iscommonshare(price_fl, shrcd);
price_fl   = price_fl(idx,:);

save('results\master.mat', 'master')
save('results\price_fl.mat','price_fl')

% Capitalizations
cap = getMktCap(master.mst.Permno, master.mst.Date,false,false,1);
cap = struct('Permnos', {getVariableNames(cap(:,2:end))}, ...
    'Dates', cap{:,1},...
    'Data', cap{:,2:end});

save('results\cap.mat','cap')

% NYSE breakpoints
bpoints = importFrenchData('ME_Breakpoints_TXT.zip','..\results');