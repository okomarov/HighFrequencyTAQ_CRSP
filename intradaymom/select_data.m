
OPT_BAD_PRICE_MULT = 10;
OPT_LAGDAY         = 1;
%% Import data

% Index data
datapath = '..\data\TAQ\sampled\5min\nobad_vw';
master   = load(fullfile(datapath,'master'),'-mat');
master   = addPermno(master.mst);
master   = master(master.Permno ~= 0,:);
master   = sortrows(master,{'Permno','Date'});

% Common shares
idx    = iscommonshare(master);
master = master(idx,:);

% Minobs 
try 
    res = loadresults('isEnoughObs');
catch
    opt = struct('Time',120000,'Minobs',30);
    res = AnalyzeImom('isEnoughObs',[],[],'data\TAQ\count\',[],[],opt);
end
[~,pos] = ismembIdDate(master.Id, master.Date, res.Id, res.Date);
idx     = res.HasEnoughObs(pos,:);
master  = master(idx,:);

% Incomplete days
idx    = isprobdate(master.Date);
master = master(~idx,:);

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

% Sample VWAP
try
    vwap = loadresults('VWAP','..\results');
catch
    mst = selectAndFilterTrades(OPT_BAD_PRICE_MULT);
    if isempty(OPT_BAD_PRICE_MULT)
        mst = mst(:, {'File','Id','Permno','Date','Isbadday'});
    else
        mst = mst(:, {'File','Id','Permno','Date','MedPrice', 'Isbadday'});
    end
    edges = [93000,100000; 120000,123000; 123000,130000; 153000,160000];
    vwap  = Analyze('VWAP',[],mst,[],[], [],struct('edgesVWAP',edges,'BadPriceMultiplier',OPT_BAD_PRICE_MULT));
end
[~,pos] = ismembIdDate(master.Permno, master.Date, vwap.Permno, vwap.Date);
vwap    = vwap(pos,:);

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
save('results\vwap.mat','vwap')
save('results\dsfquery.mat','crsp')
save('results\bpoints.mat','bpoints')

%% Half hour returns

[~,pos]           = ismembIdDate(master.Permno, master.Date, price_fl.Permno, price_fl.Date);
master.FirstPrice = price_fl.FirstPrice(pos);
master.LastPrice  = price_fl.LastPrice(pos);
master            = cache2cell(master,master.File);

ranges = [930, 1000, 1030, 1100, 1130, 1200, 1230, 1300, 1330, 1400, 1430,...
            1500,1530,1600]';
ranges = [ranges(1:end-1), ranges(2:end)]*100;
for r = 1:size(ranges,1)
    opt           = struct('HalfHourRange',ranges(r,:));
    [~, filename] = AnalyzeImom('halfHourRet',[],master,'data\TAQ\sampled\5min\nobad_vw',1,[],opt);
    
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