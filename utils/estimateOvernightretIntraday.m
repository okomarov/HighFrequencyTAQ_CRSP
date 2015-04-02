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
idx      = dsfquery.Shrcd == 10 | dsfquery.Shrcd == 11;
dsfquery = dsfquery(idx,:);
dsfquery = dsfquery(~isnan(dsfquery.Ret),{'Permno','Date','Ret'});

% % Winsorize returns at 0.1 and 99.9%
% ptiles = prctile(dsfquery.Ret,[0.1,99.9]);
% % boxplot(dsfquery.Ret)
% % idx    = in(dsfquery.Ret, ptiles);
% % boxplot(dsfquery.Ret(idx))
% dsfquery.Ret(dsfquery.Ret < ptiles(1)) = ptiles(1);
% dsfquery.Ret(dsfquery.Ret > ptiles(2)) = ptiles(2);

% Load UnID - Date pairs and keep common shares only
uniqueID = loadresults('uniqueID');
uniqueID = uniqueID(iscommonshare(uniqueID(:,{'UnID','Date'})),:);

% Map permno
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Map totret to uniqueID (by combined key to speed up)
keyA           = uint64(uniqueID.Permno) * 1e8 + uint64(uniqueID.Date);
keyB           = uint64(dsfquery.Permno) * 1e8 + uint64(dsfquery.Date);
[idx,pos]      = ismember(keyA, keyB);
uniqueID       = uniqueID(idx,:);
uniqueID.RetCC = NaN(size(uniqueID.Permno));
uniqueID.RetCC = dsfquery.Ret(pos(idx));
clear dsfquery

% Add median price
res               = loadresults('medianprice');
keyA              = uint64(uniqueID.Id) * 1e8 + uint64(uniqueID.Date);
keyB              = uint64(res.Id)      * 1e8 + uint64(res.Date);
[~,pos]           = ismember(keyA, keyB);
uniqueID.MedPrice = res.MedPrice(pos);

% Add File, From and To
path2data     = '.\data\TAQ\';
master        = load(fullfile(path2data,'master'),'-mat');
keyB          = uint64(master.mst.Id) * 1e8 + uint64(master.mst.Date);
[~,pos]       = ismember(keyA, keyB);
uniqueID.File = master.mst.File(pos);
uniqueID.From = master.mst.From(pos);
uniqueID.To   = master.mst.To(pos);
clear master
uniqueID(:,{'Permno'}) = [];

res = Analyze('return_overnight_intraday',[],uniqueID,path2data);
end