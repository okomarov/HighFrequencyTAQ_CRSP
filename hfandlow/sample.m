%% Select data
OPT_LAGDAY = 1;

% Index data
datapath = '..\data\TAQ\';
master   = load(fullfile(datapath,'master'),'-mat');
master   = addPermno(master.mst);
master   = master(master.Permno ~= 0,:);
master   = sortrows(master,{'Permno','Date'});

% Common shares
idx    = iscommonshare(master);
master = master(idx,:);

% Incomplete days
idx    = isprobdate(master.Date);
master = master(~idx,:);

% Minobs
res            = loadresults('countBadPrices','..\results');
[~,pos]        = ismembIdDate(master.Id, master.Date, res.Id, res.Date);
master.Nbadtot = res.Nbadtot(pos,:);
isEnoughObs    = (master.To-master.From+1 - master.Nbadtot) >= 79;
isEnoughObs    = [false(OPT_LAGDAY,1); isEnoughObs(1:end-OPT_LAGDAY)];
master         = master(isEnoughObs,:);

% Count
tmp       = sortrows(unstack(master(:,{'Date','Permno','File'}), 'File','Permno'),'Date');
tmp       = tmp{:,2:end};
count_all = sum(tmp~=0,2);

% CRSP returns
dsf       = loadresults('dsfquery','..\results');
[idx,pos] = ismembIdDate(dsf.Permno, dsf.Date,master.Permno, master.Date);
dsf       = dsf(idx,:);

% Beta components
beta      = loadresults('betacomponents5mon');
[idx,pos] = ismembIdDate(beta.Permno, beta.Date,master.Permno, master.Date);
beta      = beta(idx,:);

save('results\dsf.mat','dsf')
save('results\master.mat','master')
save('results\beta5minon.mat','beta')

importFrenchData('F-F_Research_Data_5_Factors_2x3_daily_TXT.zip','results');