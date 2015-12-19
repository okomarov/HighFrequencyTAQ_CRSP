% Num issues matched
path2data = '.\data\TAQ';
master    = load(fullfile(path2data, 'master'), '-mat');
master    = addPermno(master.mst);

[dates,~,subs] = unique(master.Date);
nobs           = master.To-master.From+1;
tot            = accumarray(subs,nobs);
idx            = master.Permno~=0;
matched        = accumarray(subs(idx),nobs(idx));

mean(matched./tot)
min(matched./tot)
%% 

OPT_LAGDAY = 1;

% Index data
datapath = 'data\TAQ\';
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
res            = loadresults('countBadPrices');
[~,pos]        = ismembIdDate(master.Id, master.Date, res.Id, res.Date);
master.Nbadtot = res.Nbadtot(pos,:);
isEnoughObs    = master.To-master.From+1 - master.Nbadtot >= 7;
isEnoughObs    = [false(OPT_LAGDAY,1); isEnoughObs(1:end-OPT_LAGDAY)];
master         = master(isEnoughObs,:);

% Count
tmp       = sortrows(unstack(master(:,{'Date','Permno','File'}), 'File','Permno'),'Date');
tmp       = tmp{:,2:end};
count_all = sum(tmp~=0,2);

% Filter microcaps
dsf                      = loadresults('dsfquery');
[idx,pos]                = ismembIdDate(dsf.Permno, dsf.Date,master.Permno, master.Date);
master.Price(pos(idx),1) = abs(dsf.Prc(idx));
idx                      = isMicrocap(master,'Price',OPT_LAGDAY);
master(idx,:)            = [];

% Count
tmp   = sortrows(unstack(master(:,{'Date','Permno','Price'}), 'Price','Permno'),'Date');
tmp   = tmp{:,2:end};
count = sum(~isnan(tmp),2);
plot(count)
