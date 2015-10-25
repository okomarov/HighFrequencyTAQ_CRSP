
OPT_EDGES_BAD_PRICES = [];
OPT_LAGDAY           = 1;
%% Import data

% Index data
datapath = '..\data\TAQ\sampled\5min\nobad';
master   = load(fullfile(datapath,'master'),'-mat');
master   = sortrows(master.mst,{'Permno','Date'});

% Common shares
idx    = iscommonshare(master);
master = master(idx,:);

% Has mkt cap on the previous day
cap    = getMktCap(master, OPT_LAGDAY);
idx    = cap.Cap ~= 0;
cap    = cap(idx,:);
master = master(idx,:);

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
    price_fl = Analyze('sampleFirstLast',[],mst,[],[], struct('edges',OPT_EDGES_BAD_PRICES));
end
[idx,pos] = ismembIdDate(master.Permno, master.Date, price_fl.Permno, price_fl.Date);
price_fl  = price_fl(pos(idx),:);

% Sample VWAP
try
    vwap = loadresults('VWAP','..\results');
catch
    mst = selectAndFilterTrades(OPT_EDGES_BAD_PRICES);
    if isempty(OPT_EDGES_BAD_PRICES)
        mst = mst(:, {'File','Id','Permno','Date','Isbadday','Isfewobs'});
    else
        mst = mst(:, {'File','Id','Permno','Date','MedPrice', 'Isbadday','Isfewobs'});
    end
    edges = [93000,100000; 120000,123000; 123000,130000; 153000,160000];
    vwap  = Analyze('VWAP',[],mst,[],[], struct('edgesVWAP',edges,'edges',OPT_EDGES_BAD_PRICES));
end
[idx,pos] = ismembIdDate(master.Permno, master.Date, vwap.Permno, vwap.Date);
vwap      = vwap(pos(idx),:);

% CRSP 
crsp      = loadresults('dsfquery','../results');
[idx,pos] = ismembIdDate(master.Permno, master.Date, crsp.Permno, crsp.Date);
crsp      = crsp(pos(idx),:);

save('results\master.mat', 'master')
save('results\price_fl.mat','price_fl')
save('results\vwap.mat','vwap')
save('results\dsfquery.mat','crsp')
save('results\FF49.mat','industry')