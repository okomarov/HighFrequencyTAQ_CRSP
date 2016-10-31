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
%% Count sample

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
isEnoughObs    = (master.To-master.From+1 - master.Nbadtot) >= 79;
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
%% Outlier example

% Bad cases examples
% [v,I]  = max(abs(fl.LastPrice./fl.FirstPrice - fl_bad.LastPrice./fl_bad.FirstPrice));
% tmp = getTaqData('permno',68161,19990331,19990331);
% tmp = getTaqData('permno',77010,20031030,20031031);
date = 20031031;

% Use open-to-close return example
fl    = loadresults('sampleFirstLast');
idate = ismember(fl.Date,[serial2yyyymmdd(yyyymmdd2serial(date)-1)  date]);
fl    = fl(idate,:);
fl    = addPermno(fl);

% Common and no microcaps
fl = fl(iscommonshare(fl),:);
fl = fl(~isMicrocap(fl,'LastPrice',1),:);

% Find case with bad prices
fl_bad = loadresults('sampleFirstLast_withBad');
fl_bad = fl_bad(fl_bad.Date == date,:);
fl_bad = fl_bad(ismember(fl_bad.Id,fl.Id),:);

% Average intraday return
mean(fl.LastPrice./fl.FirstPrice-1)
mean(fl_bad.LastPrice./fl_bad.FirstPrice-1)

% Bad prices for common shares only on 288 dates and count mostly < 50
bp              = loadresults('countBadPrices');
bp              = bp(bp.Nbadtot > bp.Nbadsel,:);
bp              = addPermno(bp);
bp              = bp(iscommonshare(bp),:);
bp.Nbad         = bp.Nbadtot - bp.Nbadsel;
[unDate,~,subs] = unique(bp.Date);

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(yyyymmdd2datetime(unDate), accumarray(subs,bp.Nbad), 'x')
set(gca,'TickLabelInterpreter','latex')
print('countOutliers','-depsc','-r200')
%% Counts
bad    = loadresults('countBadPrices');
con    = loadresults('consolidationcounts');
master = load('data\taq\master', '-mat');

tot      = master.mst.To - master.mst.From + 1;
obscount = [bad.Nbadsel, bad.Nbadtot-bad.Nbadsel, con.Nconsolidated, tot-bad.Nbadtot-con.Nconsolidated];

% Monthly
[unM,~,subs] = unique(con.Date/100);
tot          = accumarray(subs, tot);
obscount     = arrayfun(@(x) accumarray(subs, obscount(:,x)),1:4,'un',0);
obscount     = [obscount{:}];
prop         = bsxfun(@rdivide, obscount,tot);

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.4],'PaperPositionMode','auto')
plot(yyyymmdd2datetime(unM*100+1),tot)
set(gca,'YTick',[0,4,8,12]*1e8,'Ylim',[0,12e8],'YTickLabel',{'0','400M','800M','1200M'})
set(gca,'TickLabelInterpreter','latex')
print('obscount','-depsc','-r200','-loose')

xtick  = get(gca,'Xtick');
xtickl = get(gca,'XTickLabel');
gcapos = get(gca,'Position');

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.4],'PaperPositionMode','auto')
h = area(yyyymmdd2serial(unM*100+1),prop*100);
set(h,{'FaceColor'},...
    {[0    0.4470    0.7410]
     [0.4940    0.1840    0.5560]
     [0.9290    0.6940    0.1250]
     [0.8500    0.3250    0.0980]})    
axis tight
set(gca,'Xtick',xtick,'XTickLabel',xtickl,'Position',gcapos)
set(gca,'TickLabelInterpreter','latex')
l = legend('irregular','outliers','consolidated','good');
set(l,'interpreter','latex','location','northwest')

print('obsprop','-depsc','-r200','-loose')
%% Max num of trades per second
try
    res = loadresults('maxRecordsPerSec');
catch
    res = Analyze('maxRecordsPerSec');
end

idx              = isprobdate(res.Date) | ismember(res.Date, yyyymmdd2serial([19940404, 19961030, 19980409]));
res              = res(~idx,:);
[unDates,~,subs] = unique(res.Date);
plotdts          = serial2datetime(unDates);
maxcount         = accumarray(subs, res.Records, [],@max);

figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(plotdts, maxcount)
set(gca,'TickLabelInterpreter','latex')
print('countMaxTradeSec','-depsc','-r200')

%% Sample first and last
BAD_MULTIPLIER = 2;
mst = prepareMst();
Analyze('sampleFirstLast',[],mst,[],[],[],BAD_MULTIPLIER);
%% Sample VVWAP
BAD_MULTIPLIER = 2;
mst = prepareMst();

% Half hours
ranges = [ 930, 1000, 1030, 1100, 1130, 1200, 1230, 1300, 1330, 1400, 1430,...
          1500, 1530, 1600]'*100;
ranges = [ranges(1:end-1), ranges(2:end)];

opt             = struct('edgesVWAP',ranges,'BadPriceMultiplier',BAD_MULTIPLIER);
[out, filename] = Analyze('VWAP',[],mst,[],[],[],opt);

% Save each column separately
for r = 1:size(ranges,1)
    fname = fullfile('results', 'vwap', regexprep(filename, '.mat', sprintf('_30_%d.mat',ranges(r,1))));
    res = out(:,[1,2,r+2,end]);
    save(fname,'res')
end

% minute
ranges = serial2hhmmss(datenum(0,0,0,9,30:5:450,0)');
ranges = [ranges(1:end-1), ranges(2:end)];

for r = 1:size(ranges,1)
    opt           = struct('edgesVWAP',ranges(r,:),'BadPriceMultiplier',BAD_MULTIPLIER);
    [~, filename] = Analyze('VWAP',[],mst,[],[],[],opt);

    % Rename file
    oldName = fullfile('results',filename);
    newName = fullfile('results', 'vwap', regexprep(filename, '.mat', sprintf('_5_%d.mat',ranges(r,1))));
    movefile(oldName,newName)
end