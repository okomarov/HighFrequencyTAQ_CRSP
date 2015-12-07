%%
OPT_BAD_PRICE_MULT = 10;
OPT_LAGDAY         = 1;
%% Import data

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
isEnoughObs    = master.To-master.From+1 - master.Nbadtot >= 78;
isEnoughObs    = [false(OPT_LAGDAY,1); isEnoughObs(1:end-OPT_LAGDAY)];
master         = master(isEnoughObs,:);
% % Number of observations per 30 min
% res         = AnalyzeImom('avgObsPerBucket',[],[],'data\TAQ\count\',[],8);
% [un,~,subs] = unique(res(:,{'Date','HHMMSS'}));
% un.Avg      = accumarray(subs, res.Sum)./accumarray(subs, res.N);

% Has mkt cap on the previous day
cap    = getMktCap(master, OPT_LAGDAY);
idx    = cap.Cap ~= 0;
cap    = cap(idx,:);
master = master(idx,:);

% Sample first and last price
price_fl = loadresults('sampleFirstLast','..\results');
price_fl = addPermno(price_fl);
[~,pos]  = ismembIdDate(master.Permno, master.Date, price_fl.Permno, price_fl.Date);
price_fl = price_fl(pos,:);

% CRSP
crsp    = loadresults('dsfquery','../results');
[~,pos] = ismembIdDate(master.Permno, master.Date, crsp.Permno, crsp.Date);
crsp    = crsp(pos,:);

% NYSE breakpoints
try
    bpoints = loadresults('ME_breakpoints_TXT','..\results');
catch
    bpoints = importFrenchData('ME_Breakpoints_TXT.zip','..\results');
end
idx     = ismember(bpoints.Date, unique(cap.Date)/100);
bpoints = bpoints(idx,{'Date','Var3'});

save('results\master.mat', 'master')
save('results\price_fl.mat','price_fl')
save('results\dsfquery.mat','crsp')
save('results\bpoints.mat','bpoints')

% Check number of stocks over time
tmp = sortrows(unstack(master(:,{'Date','Permno','File'}),'File','Permno'),'Date');
plot(yyyymmdd2datetime(tmp.Date), sum(tmp{:,2:end}~=0,2));

%% Half hour returns
clearvars -except master price_fl
[~,pos]           = ismembIdDate(master.Permno, master.Date, price_fl.Permno, price_fl.Date);
master.FirstPrice = price_fl.FirstPrice(pos);
master.LastPrice  = price_fl.LastPrice(pos);
master            = cache2cell(master,master.File);

ranges = [930, 1000, 1030, 1100, 1130, 1200, 1230, 1300, 1330, 1400, 1430,...
            1500,1530,1600]';
ranges = [ranges(1:end-1), ranges(2:end)];

for r = 1:size(ranges,1)
    opt           = struct('HalfHourRange',ranges(r,:)*100);
    [~, filename] = AnalyzeImom('halfHourRet',[],master,'data\TAQ\sampled\5min\nobad_vw',[],[],opt);
    
    % Rename file
    oldName = fullfile('results',filename);
    newName = fullfile('results', regexprep(filename, '.mat', sprintf('%d.mat',ranges(r,1))));
    movefile(oldName,newName)
end

for r = 1:size(ranges,1)
    opt           = struct('HalfHourRange',ranges(r,:));
    [~, filename] = AnalyzeImom('halfHourVol',[],master,'data\TAQ\sampled\5min\nobad_vw',[],[],opt);
    
    % Rename file
    oldName = fullfile('results',filename);
    newName = fullfile('results', regexprep(filename, '.mat', sprintf('%d.mat',ranges(r,1))));
    movefile(oldName,newName)
end