
OPT_EDGES_BAD_PRICES = [];
OPT_LAGDAY           = 1;
%% Import data

% Index data
datapath = '..\data\TAQ\sampled\5min\nobad';
master   = load(fullfile(datapath,'master'),'-mat');

% Common shares
idx        = iscommonshare(master.mst);
master.mst = master.mst(idx,:);

% Has mkt cap
cap        = loadresults('mktcap','..\results');
idx        = ismembIdDate(master.mst.Permno, master.mst.Date, cap.Permno, cap.Date);
master.mst = master.mst(idx,:);

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
[idx,pos] = ismembIdDate(master.mst.Permno, master.mst.Date, price_fl.Permno, price_fl.Date);
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
[idx,pos] = ismembIdDate(master.mst.Permno, master.mst.Date, vwap.Permno, vwap.Date);
vwap      = vwap(pos(idx),:);

% CRSP 
crsp      = loadresults('dsfquery','../results');
[idx,pos] = ismembIdDate(master.mst.Permno, master.mst.Date, crsp.Permno, crsp.Date);
crsp      = crsp(pos(idx),:);

save('results\master.mat', 'master')
save('results\price_fl.mat','price_fl')
save('results\vwap.mat','vwap')
save('results\dsfquery.mat','crsp')

% Capitalizations
cap = sortrows(unstack(cap,'Cap','Permno'),'Date');
cap = struct('Permnos', {getVariableNames(cap(:,2:end))}, ...
    'Dates', cap{:,1},...
    'Data', cap{:,2:end});
% Date intersection
idx         = ismember(cap.Dates, unique(master.mst.Date));
cap.Dates   = cap.Dates(idx);
cap.Data    = cap.Data(idx,:);
% Permno expansion
xpermnos    = matlab.internal.table.numberedNames('x',unique(master.mst.Permno),false);
[idx,pos]   = ismember(cap.Permnos, xpermnos);
data        = NaN(numel(cap.Dates), numel(xpermnos));
data(:,pos(idx)) = cap.Data(:,idx);
cap.Data    = data;
cap.Permnos = xpermnos;

save('results\cap.mat','cap')

% NYSE breakpoints
bpoints = importFrenchData('ME_Breakpoints_TXT.zip','..\results');

% FF 49 industry codes
industry = getFF49IndustryCodes(master.mst);
save('results\FF49.mat','industry')