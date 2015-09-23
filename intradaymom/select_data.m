OPT_EDGES_BAD_PRICES = [];

%% Import data

% Index data
datapath = '..\data\TAQ\sampled\5min\nobad';
master   = load(fullfile(datapath,'master'),'-mat');

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


% Common shares
shrcd      = loadresults('shrcd','../results');
idx        = iscommonshare(master.mst, shrcd);
master.mst = master.mst(idx,:);
idx        = iscommonshare(price_fl, shrcd);
price_fl   = price_fl(idx,:);

save('results\master.mat', 'master')
save('results\price_fl.mat','price_fl')

% Capitalizations
cap         = getMktCap(master.mst.Permno,[],false,false,1);
cap         = struct('Permnos', {getVariableNames(cap(:,2:end))}, ...
    'Dates', cap{:,1},...
    'Data', cap{:,2:end});
% Date intersection
idx         = ismember(cap.Dates, unique(master.mst.Date));
cap.Dates   = cap.Dates(idx);
cap.Data    = cap.Data(idx,:);
% Permno expansion
xpermnos    = matlab.internal.table.numberedNames('x',unique(master.mst.Permno),false);
[~,pos]     = ismember(cap.Permnos, xpermnos);
data        = NaN(numel(cap.Dates), numel(xpermnos));
data(:,pos) = cap.Data;
cap.Data    = data;
cap.Permnos = xpermnos;

save('results\cap.mat','cap')

% NYSE breakpoints
bpoints = importFrenchData('ME_Breakpoints_TXT.zip','..\results');