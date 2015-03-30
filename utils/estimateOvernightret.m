function res = estimateOvernightret

% Load Dsfquery and filter for share type 10, 11
dsfquery = loadresults('dsfquery');
idx      = dsfquery.Shrcd == 10 | dsfquery.Shrcd == 11;
dsfquery = dsfquery(idx,:);
dsfquery = dsfquery(~isnan(dsfquery.Ret),{'Permno','Date','Ret'});

% Winsorize returns at 0.1 and 99.9%
ptiles = prctile(dsfquery.Ret,[0.1,99.9]);
% boxplot(dsfquery.Ret)
% idx    = in(dsfquery.Ret, ptiles);
% boxplot(dsfquery.Ret(idx))
dsfquery.Ret(dsfquery.Ret < ptiles(1)) = ptiles(1);
dsfquery.Ret(dsfquery.Ret > ptiles(2)) = ptiles(2);

% Load UnID - Date pairs and keep common shares only
uniqueID        = loadresults('uniqueID');
% Keep common shares only
uniqueID        = uniqueID(iscommonshare(uniqueID(:,{'UnID','Date'})),:);
% Map permno
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Map totret to uniqueID (by combined key to speed up)
keyA                 = uint64(uniqueID.Permno) * 1e8 + uint64(uniqueID.Date);
keyB                 = uint64(dsfquery.Permno) * 1e8 + uint64(dsfquery.Date);
[idx,pos]            = ismember(keyA, keyB);
uniqueID.Totret      = NaN(size(uniqueID.Permno));
uniqueID.Totret(idx) = dsfquery.Ret(pos(idx));
clear dsfquery

% Cache
path2data     = '.\data\TAQ';
master        = load(fullfile(path2data,'master'),'-mat');
keyA          = uint64(uniqueID.Id)   * 1e8 + uint64(uniqueID.Date);
keyB          = uint64(master.mst.Id) * 1e8 + uint64(master.mst.Date);
[~,pos]       = ismember(keyA, keyB);
uniqueID.File = master.mst.File(pos);
clear master
uniqueID(:,{'Permno'}) = [];

res = Analyze('return_overnight',[],uniqueID);
end