function estimateOvernightret

% Load UnID - Date pairs
uniqueID = loadresults('uniqueID');
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Load Dfequery
dfequery = loadresults('dfequery');
dfequery = dfequery(~isnan(dfequery.Ret),{'Permno','Date','Ret'});

[idx,pos] = ismemberb(uniqueID(:,{'Date','Permno'}), dfequery(:,{'Date','Permno'}),[2,3]);

end