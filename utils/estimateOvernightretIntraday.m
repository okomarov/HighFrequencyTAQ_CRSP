function res = estimateOvernightretIntraday
% Calulate overnight return from CRSP's close-to-close return and
% TAQ's open-to-close (intraday) return.

% Load overnight
try
    res = loadresults('return_intraday');
catch
    % Median price
    cached = loadresults('medianprice');
    
    % Add File
    path2data   = '.\data\TAQ\';
    master      = load(fullfile(path2data,'master'),'-mat');
    [~,pos]     = ismembIdDate(cached.Id, cached.Date, master.mst.Id, master.mst.Date);
    cached.File = master.mst.File(pos);
    clear master
    
    res = Analyze('return_intraday',[],cached,path2data);
end

% Load UnID - Date pairs
uniqueID = loadresults('uniqueID');

% Map permno
uniqueID.Permno = unid2permno(uniqueID.UnID);

% Add RetOC
[~,pos]        = ismembIdDate(uniqueID.Id, uniqueID.Date, res.Id, res.Date);
uniqueID.RetOC = res.RetOC(pos);

% Load dsfquery
dsfquery = loadresults('dsfquery');

% % Winsorize returns at 0.1 and 99.9%
% ptiles = prctile(dsfquery.Ret,[0.1,99.9]);
% % boxplot(dsfquery.Ret)
% % idx    = in(dsfquery.Ret, ptiles);
% % boxplot(dsfquery.Ret(idx))
% dsfquery.Ret(dsfquery.Ret < ptiles(1)) = ptiles(1);
% dsfquery.Ret(dsfquery.Ret > ptiles(2)) = ptiles(2);

% Add RetCC from dsfquery
[idx,pos]           = ismembIdDate(uniqueID.Permno, uniqueID.Date, dsfquery.Permno, dsfquery.Date);
uniqueID.RetCC      = NaN(size(uniqueID.Permno));
uniqueID.RetCC(idx) = dsfquery.Ret(pos(idx));
clear dsfquery

% Overnight return
uniqueID.RetCO = log((1 + uniqueID.RetCC)./(1 + uniqueID.RetOC )); 

% Filter out NaNs
idx = ~(isnan(uniqueID.RetCC) | isnan(uniqueID.RetOC));
res = uniqueID(idx,:);

filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'return_intraday_overnight');
save(fullfile('.\results\',filename), 'res')
end