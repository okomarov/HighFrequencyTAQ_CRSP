function res = estimateOvernightretIntraday
% Calulate overnight return from CRSP's close-to-close return and
% TAQ's open-to-close (intraday) return.
% NOTE: uses sampled data since it's clean.

% Load overnight
try
    res = loadresults('return_intraday');
catch
    % Median price
    cached = loadresults('medianprice');
    
    % Add File
    path2data   = '.\data\TAQ\';
    master      = load(fullfile(path2data,'master'),'-mat');
    keyA        = uint64(cached.Id) * 1e8 + uint64(cached.Date);
    keyB        = uint64(master.mst.Id) * 1e8 + uint64(master.mst.Date);
    [~,pos]     = ismember(keyA, keyB);
    cached.File = master.mst.File(pos);
    clear master
    
    res = Analyze('return_intraday',[],cached,path2data);
end

% Load UnID - Date pairs and keep common shares only
uniqueID = loadresults('uniqueID');
uniqueID = uniqueID(iscommonshare(uniqueID(:,{'UnID','Date'})),:);

% Map permno
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Add RetOC
keyA           = uint64(uniqueID.Id) * 1e8 + uint64(uniqueID.Date);
keyB           = uint64(res.Id)      * 1e8 + uint64(res.Date);
[~,pos]        = ismember(keyA, keyB);
uniqueID.RetOC = res.RetOC(pos);

% Load dsfquery (filtering for common shares done indirectly)
dsfquery = loadresults('dsfquery');

% % Winsorize returns at 0.1 and 99.9%
% ptiles = prctile(dsfquery.Ret,[0.1,99.9]);
% % boxplot(dsfquery.Ret)
% % idx    = in(dsfquery.Ret, ptiles);
% % boxplot(dsfquery.Ret(idx))
% dsfquery.Ret(dsfquery.Ret < ptiles(1)) = ptiles(1);
% dsfquery.Ret(dsfquery.Ret > ptiles(2)) = ptiles(2);

% Add RetCC from dsfqwuery
keyA                = uint64(uniqueID.Permno) * 1e8 + uint64(uniqueID.Date);
keyB                = uint64(dsfquery.Permno) * 1e8 + uint64(dsfquery.Date);
[idx,pos]           = ismember(keyA, keyB);
uniqueID.RetCC      = NaN(size(uniqueID.Permno));
uniqueID.RetCC(idx) = dsfquery.Ret(pos(idx));
clear dsfquery

% Overnight return
uniqueID.RetCO = log((1 + uniqueID.RetCC)./(1 + uniqueID.RetOC )); 

% Filter out NaNs
idx = ~(isnan(uniqueID.RetCC) | isnan(uniqueID.RetOC));
res = uniqueID(idx,:);

filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'returns_intraday_overnight');
save(fullfile('.\results\',filename), 'res')
end