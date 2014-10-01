function res = estimateOvernightret

% Load UnID - Date pairs
uniqueID        = loadresults('uniqueID');
% Keep common shares only
uniqueID        = uniqueID(iscommonshare(uniqueID(:,{'UnID','Date'})),:);
% Map permno
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Load Dfequery
dfequery = loadresults('dfequery');
% dfequery = dfequery(~isnan(dfequery.Ret),:);
dfequery = dfequery(~isnan(dfequery.Ret),{'Permno','Date','Ret'});

% Map totret to uniqueID
[idx,pos]            = ismemberb(uniqueID(:,{'Date','Permno'}), dfequery(:,{'Date','Permno'}),[2,3]);
uniqueID.Totret      = NaN(size(uniqueID.Permno));
uniqueID.Totret(idx) = dfequery.Ret(pos(idx));
% uniqueID.Open        = NaN(size(uniqueID.Permno));
% uniqueID.Close       = NaN(size(uniqueID.Permno));
% uniqueID.Open(idx)   = dfequery.Openprc(pos(idx));
% uniqueID.Close(idx)  = dfequery.Prc(pos(idx));
clear dfequery

% Cache
path2data     = '.\data\TAQ';
master        = load(fullfile(path2data,'master'),'-mat');
[~, pos]      = ismemberb(uniqueID(:,{'Date','Id'}),master.mst(:,{'Date','Id'}));
uniqueID.File = master.mst.File(pos);
clear master
uniqueID(:,{'Permno'}) = [];

res = Analyze('return_overnight',[],uniqueID);
end