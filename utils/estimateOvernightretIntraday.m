function res = estimateOvernightretIntraday
% Calulate overnight return from CRSP's close-to-close return and
% TAQ's open-to-close (intraday) return.


% Load Permno - Date pairs
mst = loadresults('masterPermno');

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
[~,pos]   = ismembIdDate(mst.Id, mst.Date, res.Id, res.Date);
mst.RetOC = res.RetOC(pos);
clear res

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
[idx,pos]      = ismembIdDate(mst.Permno, mst.Date, dsfquery.Permno, dsfquery.Date);
mst.RetCC      = NaN(size(mst.Permno));
mst.RetCC(idx) = dsfquery.Ret(pos(idx));
clear dsfquery

% Overnight return
mst.RetCO = log((1 + mst.RetCC)./(1 + mst.RetOC )); 

res      = mst;
filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'return_intraday_overnight');
save(fullfile('.\results\',filename), 'res')
end