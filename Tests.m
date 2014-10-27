%% Max trades/s
testname  = 'maxtradepsec';
try
    counts = loadresults(testname);
catch
    counts = Analyze(testname);
end
[Dates, ~, subs] = unique(counts.Date);
counts = table(yyyymmdd2datetime(Dates), accumarray(subs,counts.Maxpsec,[],@max),'VariableNames',{'Date','Maxpsec'});
plot(counts.Date, counts.Maxpsec)
ylabel 'trades/second'
print -depsc -r150 .\results\fig\maxtradepsec 
%% Counts selection rule
addpath .\utils\magnifyOnFigure\ .\utils\export_fig   

path2data = '.\data\TAQ';
testname  = 'selrulecounts';
try
    res = loadresults(testname);
catch
    res = Analyze(testname,[],[],fullfile(path2data,'T*.mat'),1);
end


refdates = serial2yyyymmdd(datenum(1993,2:209,1)-1);

% Plot 
for ii = 1:3
    data = cat(1,res{:,ii});
    % Bunch conditions
    if ii == 3
        idx           = all(data.Val == ' ' | data.Val == 'E' | data.Val == 'F' | data.Val == '@',2);
        data.Val      = cellstr(data.Val); 
        data.Val(idx) = {'  '};
    end
    data   = unstack(data,'Count','Val');
    % Filter out problematic dates
    data = data(~isprobdate(data.Date),:);
    
    % Sample dates
    data    = sampledates(data,refdates,1);
    
    vnames  = getVariableNames(data(:,2:end));
    dates   = yyyymmdd2serial(data.Date);
    data    = table2array(data(:,2:end));
    data(isnan(data)) = 0;
    [~,pos] = sort(data(1,:),'descend');
    data    = data(:,pos);
    vnames  = vnames(pos);
    
    close, figure('Color','white','Renderer', 'opengl'), colormap(lines(size(data,2)))
    % Main figure
    area(dates, bsxfun(@rdivide, data, sum(data,2))*100)
    set(gcf, 'Position', get(gcf,'Position').*[1,1,1,.4])
    dynamicDateTicks, axis tight, ylabel '%', set(gca,'Ylim',[50,100])
      
    % Legend
    switch ii
        case 1, l = {'G127 - rule compliant (0)','G127 - Display Book reported (40)'}; xMagLim = [99.9,100];
        case 2, l = 'CORR - no corrections (0)';                                       xMagLim = [97,100];
        case 3, l = {'COND - regular trades (@, ,E,F)','COND - out of sequence (Z)'};  xMagLim = [99.5,100];
    end
    legend(l,'Location','southeast');
    
    % Magnifier
    set(gca,'Units','pixels'), pos = get(gca,'Position');
    magnifyOnFigure(gca,'initialPositionMagnifier',[100,pos(2)+pos(4)-4,50,3],...
                        'initialPositionSecondaryAxes',[150,40,80,100],...
                        'secondaryAxesXLim',[728300,729000],'secondaryAxesYLim',xMagLim,...
                        'edgeWidth',1.5,'EdgeColor','white','mode','manual')
    set(gca,'XtickLabel','')

    export_fig(sprintf('.\\results\\fig\\selrulecount%d.eps',ii), '-r150','-transparent','-opengl','-a1')
end

% Plot all observations
close, figure('Color','white','Renderer', 'opengl'), colormap(lines(1))

% Main figure
data = sum(data,2);
plot(dates, data)
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,.5])
dynamicDateTicks, axis tight, ylabel 'Number of trades'
ha = gca;
ha.YRuler.Exponent = 6;
ha.YRuler.SecondaryLabel.String = 'millions';
ha.YRuler.SecondaryLabel.FontAngle = 'italic';
ha.YRuler.SecondaryLabel.Visible = 'on';

export_fig(sprintf('.\\results\\fig\\selrulecount0.eps',ii), '-r150','-transparent','-opengl','-a1')
%% Overall cleaning counts
path2data = '.\data\TAQ';
load(fullfile(path2data,'master'),'-mat')

% Map unique ID to mst
res      = loadresults('uniqueID');
[~,pos]  = ismemberb(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.UnID = res.UnID(pos);

% Bad prices days
res          = loadresults('badprices');
[~,pos]      = ismemberb(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.Nbadsel  = res.Nbadsel(pos);
mst.Nbadtot  = res.Nbadtot(pos);
mst.Isbadday = res.Isbadday(pos);

% Bad series
nobs            = mst.To - mst.From +1;
totbad          = accumarray(mst.UnID(mst.Isbadday), nobs(mst.Isbadday),[max(mst.UnID),1]);
totobs          = accumarray(mst.UnID, nobs);
badseries       = totbad./totobs > .1;
badseries(end)  = true; % for the unmatched
mst.Isbadseries = badseries(mst.UnID);

% Count losing obs with timestamp consolidation
res               = loadresults('consolidationcounts');
[~,pos]           = ismemberb(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
mst.Nconsolidated = res.Nconsolidated(pos);

% Select on basis of minimum number of observations
ngoodtrades    = mst.To-mst.From+1 - mst.Nbadtot - mst.Nconsolidated;
mst.Ifewtrades = ngoodtrades < 13;
perfew         = accumarray(mst.UnID, mst.Ifewtrades)./accumarray(mst.UnID, 1) > .5;
mst.Ifewseries = perfew(mst.UnID);

% Accumulate monthly
[undates,~,subs] = unique(mst.Date/100);
sz               = [numel(undates),1];
Nbadsel          = accumarray(subs, mst.Nbadsel);
Nbad             = accumarray(subs, mst.Nbadtot) - Nbadsel;
idx              = mst.Isbadday;
Nbadday          = accumarray(subs(idx), nobs(idx)-mst.Nbadtot(idx),sz);
idx              = mst.Isbadseries & ~mst.Isbadday;
Nbadseries       = accumarray(subs(idx), nobs(idx)-mst.Nbadtot(idx),sz);
Nconsolidated    = accumarray(subs, mst.Nconsolidated);
idx              = mst.Ifewtrades & ~(mst.Isbadday | mst.Isbadseries);
Nfewobs          = accumarray(subs(idx), nobs(idx)-mst.Nbadtot(idx)-mst.Nconsolidated(idx),sz);
idx              = mst.Ifewseries & ~ (mst.Ifewtrades  | mst.Isbadday | mst.Isbadseries);
Nfewdays         = accumarray(subs(idx), nobs(idx)-mst.Nbadtot(idx)-mst.Nconsolidated(idx),sz);
Ntot             = accumarray(subs,nobs);

% plot dates
if undates(1) < 19921231
    undates = double(undates);
    plotdates = datenum(fix(undates/100), mod(undates,100)+1, 1)-1;
else
    plotdates = yyyymmdd2serial(undates);
end

% Plot cleaning
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,.5],'PaperPositionMode','auto')

relval = bsxfun(@rdivide,[Nbad, Nbadday, Nbadseries],Ntot)*100;
colormap(lines(size(relval,2)))
area(plotdates, relval,'LineStyle','none')

dynamicDateTicks
axis tight, set(gca,'Layer','top','Ylim',[0,1]),ylabel '%'
legend({'Bad observations','Bad days','Bad series'},'Location','northwest')

print .\results\fig\cleaningcounts.eps -depsc -r300 -opengl

% Plot cleaning
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,.5],'PaperPositionMode','auto')

relval = bsxfun(@rdivide,[Nbadsel, Nbad+Nbadday+Nbadseries, Nconsolidated, Nfewobs, Nfewdays],Ntot)*100;
colormap(lines(size(relval,2)))
area(plotdates, relval,'LineStyle','none')
% mean(relval)

dynamicDateTicks
axis tight, set(gca,'Layer','top','Ylim',[0,100]),ylabel '%'
legend({'Selection','Cleaning','Consolidation','Minimum obs.','Min days'},'Location','northwest')

print .\results\fig\allrulescounts.eps -depsc -r300 -opengl
%% Display Book (G127 - 40) Keep? [YES]
path2data = '.\data\TAQ';
master    = load(fullfile(path2data, 'master'), '-mat');
s         = load(fullfile(path2data, 'T0300.mat'));

% Count DB trades by Id
[~,pmst]    = histc(find(s.data.G127_Correction(:,1) == 40), [s.mst.From; s.mst.To(end)]);
mst         = s.mst(unique(pmst),:);
[un,~,subs] = unique(s.mst.Id(pmst));
dbtrades    = table(un, accumarray(subs,1),'VariableNames', {'Id','Count'});
[~,imax]    = max(dbtrades.Count);
% imax      = 100;
disp(dbtrades(imax,:))
id          = dbtrades.Id(imax);
symbol      = s.ids{id}
date        = mst.Date(find(mst.Id == id,1,'first'));

% Plot
sample = getTaqData(master, symbol, date,date);
i40    = sample.G127_Correction(:,1) == 40;
igood  = sample.G127_Correction(:,1) == 0;
x      = hhmmssmat2serial(sample.Time);
plot(x(i40), sample.Price(i40),'xr', x(igood), sample.Price(igood),'.b',...
     x, sample.Price, '-g')
%% Counts by type
try
    TAQmaster = loadresults('TAQmaster');
catch
    TAQmaster = importMst('.\data\TAQ\raw');
end

TAQmaster = unique(TAQmaster(:, {'SYMBOL','FDATE','TYPE'}));

% Time consolidation
idx       = isfeatchange(TAQmaster(:,[1,3,2]),2:3);
TAQmaster = TAQmaster(idx,:);

% Link to number of records
master = load(fullfile(path2data, 'master'), '-mat');

% Subscripts by date
[unDates, ~, subs] = unique(master.mst.Date/100);
master.mst.Subs    = uint8(subs);

% Preallocation
master.mst.Val     = zeros(numel(subs),1,'uint8');
master.ids = regexprep(master.ids,'p','PR');

% LOOP by symbol
for ii = 1:numel(master.ids)
    symbol   = master.ids{ii};
    iTAQ     = strcmpi(TAQmaster.SYMBOL,symbol);
    if ~any(iTAQ)
        fprintf('No match: %s - iter %d\n', symbol,ii), continue
    end
    types = TAQmaster.TYPE(iTAQ);
    if types == 0, continue, end
    % Map mst records to date bins from TAQmaster
    imst     = master.mst.Id == ii;
    if numel(types) == 1
        master.mst.Val(imst) = types;
        continue
    end
    dates    = master.mst.Date(imst,:);
    datebins = [TAQmaster.FDATE(iTAQ); inf];
    [~,dmap] = histc(dates, datebins);
    if dmap(1) == 0
        fprintf('Data before type: %s - iter %d\n', symbol,ii)
        dmap(dmap == 0) = ones(1,'uint8');
    end
    % Count by month
    master.mst.Val(imst) = types(dmap);
end
% Group by month and type and count
[typecounts, ~, subs] = unique(master.mst(:,{'Subs','Val'}));
typecounts.Subs       = unDates(typecounts.Subs);
typecounts.Counts     = accumarray(subs, master.mst.To-master.mst.From+1);

% Unstack
typecounts = dataset2table(unstack(typecounts,'Counts','Val'));
vnames     = {'Common', 'Preferred', 'Warrant', 'Right', 'Other', 'Derivative'};
typecounts = setVariableNames(typecounts, ['Dates' vnames]);

% Save
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'typecounts')), 'typecounts')

% Plot
dates = datenum(double(typecounts.Dates/100), double(rem(typecounts.Dates,100)+1), 1)-1;
data  = table2array(typecounts(:,2:end));
data(isnan(data)) = 0;
subplot(211)
area(dates, data)
dynamicDateTicks
axis tight
title('Montly counts of price observations by type (absolute and %)')
l = legend({'common', 'preferred', 'warrant', 'right', 'other', 'derivative'});
set(l,'Location','NorthWest')
subplot(212)
area(dates, bsxfun(@rdivide, data, sum(data,2))*100)
dynamicDateTicks
axis tight
%% TAQ master vs data
% Check for how many symbols in master there is no actual data

% TAQmaster symbol and min date
TAQmaster = loadresults('TAQmaster');
[unsymb,~,subs] = unique(TAQmaster(:,'SYMBOL'));
unsymb.Date     = accumarray(subs, TAQmaster.FDATE,[],@min);

% Master
path2data  = '.\data\TAQ';
master     = load(fullfile(path2data, 'master'), '-mat');

[~,ia,ib] = intersect(unsymb.SYMBOL, master.ids);

% Traded but not in master records
anomalous = setdiff(1:numel(master.ids),ib); 
master.ids(anomalous)

% Did not trade but in master records
nontraded = setdiff(1:size(unsymb,1),ia); 
unsymb.SYMBOL(nontraded)
%% CRSP link coverage

% TAQ2CRSP link
taq2crsp = loadresults('taq2crsp');
taq2crsp = taq2crsp(~isnan(taq2crsp.permno),:);
taq2crsp = sortrows(taq2crsp,{'symbol','datef'});

% Master file
path2data  = '.\data\TAQ';
master     = load(fullfile(path2data, 'master'), '-mat');

% Subscripts by date
[unDates, ~, subs] = unique(master.mst.Date/100);
master.mst.Subs    = uint8(subs);

% Preallocation
master.ids       = regexprep(master.ids,'p','PR');
master.ids       = regexprep(master.ids,'\.','');
master.mst.Score = zeros(numel(subs),1,'uint8');

% 
[unsymbols,~,subsymbols] = unique(taq2crsp.symbol);
for ii = 1:numel(unsymbols)
    % taq2crsp info
    symbol = unsymbols{ii};
    itaq   = subsymbols == ii;
    dates  = [taq2crsp.datef(itaq);inf];
    permnos = uint8([0; taq2crsp.score(itaq)]);
    
    % master info
    id = find(strcmpi(master.ids,symbol));
    if isempty(id)
        fprintf('No match: %s - iter %d\n', symbol,ii), continue
    end
    if numel(id) == 1
        imst = master.mst.Id == id;
    else
        imst = ismember(master.mst.Id,id);
    end
    tmp  = master.mst.Date(imst);
    [~, datebin] = histc(tmp, dates);
    master.mst.Score(imst) = permnos(datebin+1);
end

% Group by month and score and count
[counts, ~, subs] = unique(master.mst(:,{'Subs','Score'}));
counts.Subs       = unDates(counts.Subs);
counts.Counts     = accumarray(subs, master.mst.To-master.mst.From+1);

% Unstack
counts   = dataset2table(unstack(counts,'Counts','Score'));
vnames = {'0. Unmatched', '1. CUSIP', '2. 1 + expand exact name',...
          '3. 2 + expand through new cusip', '4. SYMBOL + DATE',...
          '5. SYMBOL + lev(NAME)','6. 5 + expand through new cusip',...
          '7. lev(NAME)', '8. 7 + expand through new cusip'};
counts.Properties.VariableNames{1} = 'Dates';

% Save
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'matchcounts')), 'counts')

% Plot
dates = datenum(double(counts.Dates/100), double(rem(counts.Dates,100)+1), 1)-1;
data  = table2array(counts(:,2:end));
data(isnan(data)) = 0;
subplot(211)
area(dates, data)
dynamicDateTicks
axis tight
title('Montly counts of price observations by taq2crsp match (absolute and %)')
l = legend(vnames);
set(l,'Location','NorthWest')
subplot(212)
area(dates, bsxfun(@rdivide, data, sum(data,2))*100)
dynamicDateTicks
axis tight
%% SHRCD selection/counts

% Load msenames
try
    res = loadresults('shrcd');
catch
    res = mapShrcd2mst;
end
% Master file
path2data  = '.\data\TAQ';
load(fullfile(path2data, 'master'), '-mat');

% Counts 
% -------------------------------------------------------------------------
% Subscripts by date
[unDates, ~, Subs] = unique(mst.Date/100);
res.Subs = uint8(Subs);

% Group by month and score and count

[counts, ~, subs] = unique(res(:,{'Subs','Shrcd'}));
counts.Subs       = unDates(counts.Subs);
counts.Counts     = accumarray(subs, mst.To-mst.From+1);

% Unstack
counts   = unstack(counts,'Counts','Shrcd');
counts.Properties.VariableNames{1} = 'Dates';
vnames = getVariableNames(counts);

% Select a few
map = table({'not matched';'common-undefined';'common';'common-incorporated not US';'ADR';'ETFs'},'RowNames',{'x0','x10','x11','x12','x31','x73'});
idx = ismember(vnames(2:end),  map.Properties.RowNames);

% Plot
dates = datenum(double(counts.Dates/100), double(rem(counts.Dates,100)+1), 1)-1;
data  = table2array(counts(:,2:end));
data  = [data(:,idx), nansum(data(:,~idx),2)];
data(isnan(data)) = 0;
subplot(211)
area(dates, data)
dynamicDateTicks
axis tight
title('Montly counts of price observations by Share Type Code (absolute and %)')
l = legend([map.Var1; 'other']);
set(l,'Location','NorthWest')
subplot(212)
area(dates, bsxfun(@rdivide, data, sum(data,2))*100)
dynamicDateTicks
axis tight

% Save
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'shrcdcounts')), 'counts')

%% Count 0 returns

testname = 'countnullrets';
try
    counts = loadresults(testname);
catch
    path2data = '.\data\TAQ\sampled';
    counts    = Analyze(testname,[],[], fullfile(path2data,'S5m_*.mat'),1);
end

% All
[refdates,~,subs] = unique(counts.Date);
avgcounts         = accumarray(subs, counts.Nullrets,[],@mean);

subplot(311)
plotdates = datetime(yyyymmdd2serial(refdates),'ConvertFrom','datenum');
plot(plotdates,avgcounts)
title 'ALL: average # of null returns out of 78 daily returns'

% Only sp500 members
isp500    = issp500member(counts(:,{'Date','UnID'}));
subs      = ismembc2(counts.Date(isp500),refdates);
avgcounts = accumarray(subs, counts.Nullrets(isp500),[],@mean);

subplot(312)
plot(plotdates,avgcounts)
title 'SP500 members'

% Spyders
spy         = sortrows(counts(counts.UnID == 29904,{'Date','Nullrets'}),'Date');

counts      = NaN(numel(refdates),1);
pos         = ismembc2(spy.Date,refdates);
counts(pos) = spy.Nullrets;

subplot(313)
plot(plotdates,counts)
title 'SPY'

saveas(gcf,'.\results\NullRetCounts.png')
%% Count 0 returns on mkt cap
addpath .\utils\

% Filter settings
sp500only  = true;
commononly = true;

testname = 'countnullrets';
try
    counts = loadresults(testname);
catch
    path2data = '.\data\TAQ\sampled';
    counts    = Analyze(testname,[],[], fullfile(path2data,'S5m_*.mat'),1);
end

if commononly
    idx    = iscommonshare(counts(:,{'UnID','Date'}));
    counts = counts(idx,:);
end

if sp500only
    idx    = issp500member(counts(:,{'UnID','Date'}));
    counts = counts(idx,:);
end

% Unstack counts
unids  = unique(counts.UnID);
counts = unstack(counts(:,{'UnID','Date','Nullrets'}),'Nullrets','UnID');

% Filter out problematic dates
counts = counts(~isprobdate(counts.Date),:);

% Sample dates
refdates = serial2yyyymmdd(datenum(1993,2:209,1)-1);
counts   = sampledates(counts,refdates,1);

% Get mkt capitalizations
mktcap = estimateMktCap(unids, refdates);

% Intersect/extract data
[~,icount,icap] = intersect(getVariableNames(counts), getVariableNames(mktcap));
counts          = table2array(counts(:, icount(2:end)));
mktcap          = table2array(mktcap(:, icap(2:end)));

% Intersect NaNs
inan = isnan(mktcap) | isnan(counts);
mktcap(inan) = NaN;
counts(inan) = NaN;

% Index Betas by market cap
ptiles = prctile(mktcap,10:10:90,2);
N      = size(ptiles,1);
ptiles = [zeros(N,1), ptiles, inf(N,1)];
subs = mktcap;
for r = 1:N
    [~, subs(r,:)] = histc(mktcap(r,:), ptiles(r,:));
end
rsubs = repmat((1:N)',1,size(mktcap,2));
averages = accumarray([rsubs(~inan),subs(~inan)], counts(~inan),[N, 10],@mean);

plotdates = datetime(yyyymmdd2serial(refdates),'ConvertFrom','datenum');
plot(plotdates, averages)
legend(arrayfun(@(x) sprintf('%d^{th} ',x),10:10:100,'un',0))

if sp500only
    title 'Null ret counts: mkt cap percentile averages of SP500 members'
    saveas(gcf, '.\results\NullRetCounts_SP500member_mktcap.png')
else
    title 'Null ret counts: mkt cap percentile averages of all stocks'
    saveas(gcf, '.\results\NullRetCounts_All_mktcap.png')
end
%% Beta quantiles
betaPercentiles(10:10:90,1,5,true,false,true,true)
%% Beta quantiles on mkt cap
addpath .\utils\

% Filter settings
sp500only  = true;
commononly = true;

% Get Betas
[Betas, unids] = getBetas(sp500only, commononly);

% Filter out problematic dates
Betas = Betas(~isprobdate(Betas.Date),:);

% Sample/expand
refdates = serial2yyyymmdd(datenum(1993,2:209,1)-1);
Betas    = sampledates(Betas,refdates,1);

% Get mkt capitalizations
mktcap = estimateMktCap(unids, refdates);

% Intersect/extract data
[~,ibetas,icap] = intersect(getVariableNames(Betas), getVariableNames(mktcap));
Betas           = table2array(Betas(:, ibetas(2:end)));
mktcap          = table2array(mktcap(:, icap(2:end)));

% Intersect NaNs
inan = isnan(mktcap) | isnan(Betas);
mktcap(inan) = NaN;
Betas(inan) = NaN;

% Index Betas by market cap
ptiles = prctile(mktcap,10:10:90,2);
N      = size(ptiles,1);
ptiles = [zeros(N,1), ptiles, inf(N,1)];
subs = mktcap;
for r = 1:N
    [~, subs(r,:)] = histc(mktcap(r,:), ptiles(r,:));
end
rsubs = repmat((1:N)',1,size(Betas,2));
averages = accumarray([rsubs(~inan),subs(~inan)], Betas(~inan),[N, 10],@mean);

plotdates = datetime(yyyymmdd2serial(refdates),'ConvertFrom','datenum');
plot(plotdates, averages)
legend(arrayfun(@(x) sprintf('%d^{th} ',x),10:10:100,'un',0))

if sp500only
    title 'Betas: mkt cap percentile averages of SP500 members'
    saveas(gcf, '.\results\BetasSP500_mktcap_proxy.png')
else
    title 'Betas: mkt cap percentile averages of all stocks'
    saveas(gcf, '.\results\BetasAll_mktcap_proxy.png')
end
%% Cap-weighted Beta

% Filter settings
sp500only  = true;
commononly = false;

% Get Betas
[Betas, unids] = getBetas(sp500only, commononly);

% Sample/expand
refdates = serial2yyyymmdd(datenum(1993,2:234,1)-1);
Betas    = sampledates(Betas,refdates,1);

% Get mkt capitalizations
mktcap = estimateMktCap(unids, refdates);

% Intersect/extract data
[~,ibetas,icap] = intersect(getVariableNames(Betas), getVariableNames(mktcap));
Betas           = table2array(Betas(:, ibetas(2:end)));
mktcap          = table2array(mktcap(:, icap(2:end)));

% Weights
weights             = bsxfun(@rdivide, mktcap, nansum(mktcap,2));
Betas               = Betas.*weights;
Betas(weights == 0) = NaN;

% Plot
plotdates = datetime(yyyymmdd2serial(refdates),'ConvertFrom','datenum');
plot(plotdates, nansum(Betas,2))

if sp500only
    title 'Cap-weighted sum of SP500 Betas (with proxy)'
    saveas(gcf, '.\results\BetaCapWeightSP500_proxy.png')
else
    title 'Cap-weighted sum of all Betas (with proxy)'
    saveas(gcf, '.\results\BetaCapWeightAll_proxy.png')
end
%% SP500 momentum
addpath .\utils\
spdir  = '.\data\CRSP';

% Filter settings
sp500only  = true;
commononly = true;

% Betas spy
Betasspy = getBetas(sp500only, commononly,'20140704_1722_betas');

% Load SPY (etf)
loadname = 'spysampled';
try
    spy = loadresults(loadname);
catch
    master    = load(fullfile(path2data,'master'),'-mat');
    spy       = getTaqData(master, 'SPY',[],[],'Price',path2data);
    save(fullfile(resdir,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),loadname)), 'spy')
end

% Betasspy*spy
spy             = iprice2dret(spy);
[~,pos]         = ismember(Betasspy.Date,spy.Date);
ret             = [NaN; spy.Ret];
Betasspy.Sysret = ret(pos+1).*Betasspy.Beta;

% Get Betas
[Betas, unids] = getBetas(sp500only, commononly);

% Load sp500proxy
spproxy = loadresults('sp500proxy');

% Betas*proxy
spproxy      = iprice2dret(spproxy);
[~,pos]      = ismember(Betas.Date,spproxy.Date);
ret          = [NaN; spproxy.Ret];
Betas.Sysret = ret(pos+1).*Betas.Beta;

% Daily rets [Use the dsfquery from crsp]
try
    rest = loadresults('dailyret');
catch
    path2data = '.\data\TAQ\sampled';
    rets = Analyze('dailyret', [], [], fullfile(path2data,'S5m*.mat'));
end
[rets.Netprx,rets.Netspy] = deal(NaN(size(rets,1),1));

% Net rets proxy
[~,pos]          = ismember(Betas(:,{'UnID','Date'}), rets(:,{'UnID','Date'}));
rets.Netprx(pos) = rets.Dret(pos) - Betas.Sysret;

% Net rets spy
[~,pos]          = ismember(Betasspy(:,{'UnID','Date'}), rets(:,{'UnID','Date'}));
rets.Netspy(pos) = rets.Dret(pos) - Betasspy.Sysret;

momstrat(setVariableNames(rets(~isnan(rets.Netprx),1:4),{'UnID','Date','Dayret','Netret'}))
momstrat(setVariableNames(rets(~isnan(rets.Netspy),[1:3, 5]),{'UnID','Date','Dayret','Netret'}))

legend('R - beta*mkt','Whole returns')
title('Long/short 90th/10th percentile based on previous month performance (no overnight)')

%% Check Betas
addpath .\utils\ .\utils\nth_element\ .\utils\MFE

symbol = 'AAPL';
dovernight = false;

% Load SP500
SPY = loadresults('spysampled');

% Which files to loop for
d      = '.\data\TAQ';
load(fullfile(d, 'master'), '-mat');
idx    = mst.Id == find(strcmpi(ids,symbol));
mst    = mst(idx,:);
files  = unique(mst.File);
nfiles = numel(files);

% Sampling grid
step  = 5/(60*24);
grid  = (9.5/24:step:16/24)';
betas = cell(nfiles,1);

p = poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'}, 'debug',false);

tic
% LOOP by files
parfor (f = 1:nfiles, 4)
    
    disp(f)
    % Load data
    s = load(fullfile(d, sprintf('T%04d.mat',files(f))));
    % Select symbol
    idx   = s.mst.Id == find(strcmpi(s.ids,symbol));
    s.mst = s.mst(idx,:);
    nmst  = size(s.mst,1);
    
    res = NaN(nmst,2);
    
    % LOOP by days
    for ii = 1:nmst
        % Data selection
        from = s.mst.From(ii);
        to   = s.mst.To(ii);
        data = s.data(from:to,:);
        data = data(~selecttrades(data),{'Time','Price'});
        % Datetime
        day           = yyyymmdd2serial(s.mst.Date(ii));
        data.Datetime = day + hhmmssmat2serial(data.Time);
        % Median
        [dates, ~,subs] = unique(data.Datetime);
        prices          = accumarray(subs,data.Price,[],@fast_median);
        
        % SP500
        iday = fix(SPY.Datetime) == day;
        if ~any(iday), continue, end
        
        SPprices  = SPY.Price(iday);
        SPdates   = SPY.Datetime(iday);
        notnan    = ~isnan(SPprices);
        [~, rcss] = realized_covariance(SPprices(notnan), SPdates(notnan), prices, dates,...
            'unit','fixed', grid+day,1);
        % Overnight
        if ii ~= 1 && dovernight
            if prices(1)/s.data.Price(s.mst.To(ii-1))-1 < 0.1
                onp = prices(1)/s.data.Price(s.mst.To(ii-1))-1;
            else
                onp = 0;
            end
            pos  = find(iday,1,'first');
            onsp = SPY.Price(pos)./ SPY.Price(pos-1)-1;
        else
            onp = 0;
            onsp = 0;
        end
        res(ii,:) = [day (rcss(2)+onp*onsp)/(rcss(1)+onsp^2)];
    end
    betas{f,1} = res;
end
betas = cat(1, betas{:});
betas = betas(~isnan(betas(:,2)),:);
save debugstate2

refdates = datenum(1993,2:234,1)-1;
out = interp1(betas(:,1), betas(:,2), refdates');
plot(refdates, out)
dynamicDateTicks
hold on

% Compare against saved betas [Different!]
taq2crsp = loadresults('taq2crsp');
betas2 = loadresults('betas');
ID = taq2crsp.ID(strcmpi(taq2crsp.symbol,symbol),:);
betas2 = betas2(betas2.UnID == ID,:);
betas2 = betas2(~isnan(betas2.Beta),:);
refdates = serial2yyyymmdd(refdates);
out = interp1(double(betas2.Date), betas2.Beta, refdates');
plot(yyyymmdd2serial(refdates), out,'r')
%% Detailed check
date = 19970806;
step = 5/(60*24);
grid = (9.5/24:step:16/24)';

% Sampled
path2data = '.\data\TAQ\sampled';
master    = load(fullfile(path2data, 'master'), '-mat');
spys      = getTaqData(struct(master), 'spy', date, date,[], path2data);
aapls     = getTaqData(struct(master), 'aapl', date, date,[], path2data);

% Full
path2data = '.\data\TAQ\';
master    = load(fullfile(path2data, 'master'), '-mat');
spy       = getTaqData(struct(master), 'spy', date, date,[], path2data);
aapl      = getTaqData(struct(master), 'aapl', date, date,[], path2data);

% Select
spy  = spy(~selecttrades(spy),{'Datetime','Price'});
aapl = aapl(~selecttrades(aapl),{'Datetime','Price'});

% Median
[dates, ~,subs] = unique(spy.Datetime);
prices          = accumarray(subs,spy.Price,[],@fast_median);
spy             = table(dates,prices,'VariableNames', getVariableNames(spy));

[dates, ~,subs] = unique(aapl.Datetime);
prices          = accumarray(subs,aapl.Price,[],@fast_median);
aapl            = table(dates,prices,'VariableNames', getVariableNames(aapl));

[~, rcss] = realized_covariance(spy.Price, spy.Datetime, aapl.Price,aapl.Datetime,...
    'unit','fixed', grid+yyyymmdd2serial(date),1);
save debugstate

%% SPX (tickwrte) vs SP500 (CRSP)
d     = '.\data\';
opts  = {'Delimiter',',','ReadVarNames',1};
start = datenum('31/12/1992','dd/mm/yyyy');

% Load SP500 CRSP index
SP500crsp = dataset('File',fullfile(d,'SP500','dsp500.csv'),opts{:});
SP500crsp = replacedata(SP500crsp,@yyyymmdd2serial,'caldt');
% Keep those with enddate > 31/12/1992 and select 'caldt' and 'spindx'
SP500crsp = SP500crsp(SP500crsp.caldt > start,:);
SP500crsp = SP500crsp(:,{'caldt','spindx'});

% Load SPX (tickwrite)
filename  = unzip(fullfile(d,'Tickwrite','SP.zip'), fullfile(d,'Tickwrite'));
SP500tick = dataset('File',filename{:},opts{:},'format','%f%f%f%*[^\n]');
SP500tick = replacedata(SP500tick,@(x) yyyymmdd2serial(x + ((x>9e5).*19e6 + (x<=9e5).*20e6)),'Date');
SP500tick = replacedata(SP500tick,@(x) hhmmssfff2serial(x),'Time');
% Select close prices
SP500tick = SP500tick(diff(SP500tick.Date) > 0,:);

% Intersect dates and plot
[dates,icrsp,itick] = intersect(SP500crsp.caldt, SP500tick.Date);

% Indices
figure('pos',[300,500,800,400])
axes('pos',[.05,.12,.58,.8])
plot(dates, SP500crsp.spindx(icrsp), dates, SP500tick.Price(itick))
datetick('x')
axis tight
legend('SP500 crsp','SPX tickwrite'); legend('location','NorthWest')

% Percentage absolute tracking error (perATE)
axes
perATE = abs(tick2ret(SP500crsp.spindx(icrsp)) - tick2ret(SP500tick.Price(itick)))*100;
boxplot(gca,perATE)
set(gca,'pos',[.7,.12,.25,.65])
str = {'PERcentage Absolute Tracking Error'
    sprintf('%-10s%5.2f%%','mean:',mean(perATE))
    sprintf('%-10s%5.2f%%','std:',std(perATE))};
h = title(str,'hor','left');
oldpos = get(h,'pos');
set(h,'pos',[0.5,oldpos(2:end)])
% Save figure
print(gcf,fullfile(cd,'results','SP500 vs SPX - tracking error.png'),'-dpng','-r300')
% Clean
delete(filename{:})
%% TAQspy vs TICKWRITEspy

%% Discard first years [No]
cd C:\HFbetas
addpath .\utils\
d   = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')
load .\results\20131121_1703_avgtimestep.mat

% Number of few trades per month
ifewtrades = isnan(res.Timestep) | res.Timestep > 1/48;
[unym,~,subs] = unique(mst.Date(ifewtrades)/100);
dates = yyyymmdd2serial( double(unym)*100+1);
area(dates, accumarray(subs,1))
dynamicDateTicks
axis tight
% Concentration of few trades
plot(dates, accumarray(subs,single(mst.Id(ifewtrades)),[],@ginicoeff))
dynamicDateTicks
axis tight
% Histogram of percentage of few trades per security
perfew = accumarray(mst.Id, ifewtrades)./accumarray(mst.Id, 1);
hist(perfew,100);
%% One symbol per csv
cd C:\HFbetas
addpath .\utils .\utils\mcolon
writeto = 'E:\HFbetas\data\TAQ\Kasra';

% Open matlabpool
% if matlabpool('size') == 0
%     matlabpool('open', 4, 'AttachedFiles',{'.\utils\poolStartup.m'})
% end

setupemail
try
    tic
    d   = '.\data\TAQ';
    dd  = dir(fullfile(d,'*.mat'));
    N   = numel(dd);
    
    % Series to zip after they won't appear in the next files
    load(fullfile(d,'master'),'-mat','ids','mst')
    mst        = mst(:,{'Id','File'});
    [idx, pos] = unique(mst.Id,'last');
    tozip      = arrayfun(@(x) ids(idx(mst.File(pos)==x)),1:max(mst.File),'un',0);
    
    
    % Create all files, then append
    filenames = cellfun(@(x) sprintf('%s_.csv',x),ids,'un',0);
    existing  = dir(fullfile(writeto,'*.csv'));
    filenames = setdiff(filenames,{existing.name});
    try
        for ii = 1:numel(filenames) % 6617 24510 doesn't write!
            fid = fopen(fullfile(writeto,filenames{ii}),'w');
            fprintf(fid,'SYMBOL,DATE,TIME,PRICE,SIZE,G127,CORR,COND,EX\n');
            fclose(fid);
        end
    catch
    end
    % Update non written
    existing  = dir(fullfile(writeto,'*.csv'));
    filenames = setdiff(filenames,{existing.name});
    
    clear mst ids existing
    
    % Load file
    for f = 1:N
        disp(f)
        % Load data and slice it by Id
        s      = load(fullfile(d,dd(f).name));
        nids   = numel(s.ids);
        s.mst  = arrayfun(@(x) s.mst(s.mst.Id == x,{'Date','From','To'}),1:nids,'un',0);
        
        % LOOP by id
        for t = 1:nids
            symb = s.ids{t};
            % Avoid some files like BSp etc...
            if any(strcmp(sprintf('%s_.csv',symb),filenames))
                continue
            end
            % Open in append mode
            file2append = fullfile(writeto, sprintf('%s_.csv', symb));
            fid = fopen(file2append,'a');
            for ii = 1:size(s.mst{t},1)
                fmt = sprintf('%s,%d,%s',symb, s.mst{t}.Date(ii),'%d:%d:%d,%f,%d,%d,%d,%c%c,%c\n');
                fprintf(fid, fmt, double(s.data(s.mst{t}.From(ii):s.mst{t}.To(ii),:))');
            end
            fclose(fid);
            % Zip if no more data in next files
            if any(strcmp(symb,tozip{f}))
                % Try 3 times to zip
                for ii = 1:3
                    zipname = fullfile(writeto, sprintf('%s_.zip', symb));
                    zip(zipname,file2append)
                    try
                        valid = org.apache.tools.zip.ZipFile(java.io.File(zipname));
                    catch err
                        if ii == 3, disp(err.message), end
                    end
                end
                delete(file2append)
            end
        end
    end
    
    % Collect all results and convert to dataset
    % Export results and notify
    message = sprintf('Task ''%s'' terminated in %s','Ticker2csv',sec2time(toc)); disp(message)
    sendmail('o.komarov11@imperial.ac.uk', message,'')
catch err
    filename = fullfile(writeto, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''','Ticker2csv'), err.message, {filename})
    rethrow(err)
end
% matlabpool close
rmpref('Internet','SMTP_Password')
%% Smooth Betas [Not using]
addpath .\utils\

% Load Betas
Betas = loadresults('betas');

% Sort betas
Betas = sortrows(Betas,{'UnID','Date'});

% Index by ID
[unID, ~,subsID] = unique(Betas.UnID);

% Exclude series with less than 20 observations, i.e. one month of data
ifewdays = accumarray(subsID,1) < 20;
ikeep    = ~ismember(Betas.UnID, unID(ifewdays));

% Moving averages
% tmp = accumarray(subsID(ikeep), Betas.Date(ikeep),[],@issorted, true);
sz     = size(unID);
Betasd = dataset();
tmp    = accumarray(subsID(ikeep), Betas.UnID(ikeep),sz,@(x) {x(5:end)});
Betasd.ID = cat(1,tmp{:});
tmp    = accumarray(subsID(ikeep), Betas.Date(ikeep),sz,@(x) {x(5:end)});
Betasd.Date = cat(1,tmp{:});
tmp = accumarray(subsID(ikeep), Betas.Beta(ikeep),sz,@(x) {conv(x,ones(1,5)/5,'valid')});
Betasd.SMA = cat(1,tmp{:});
% tmp = accumarray(subsID(ikeep), Betas.Beta(ikeep),sz,@(x) {movavg(x,1,5,'e')});
% Betasd.EMA = cat(1,tmp{:});

% Weekly betas
% Remove overlapping
[un,~,subs] = unique(Betas(:,{'UnID','Date'}));
overlap     = un(accumarray(subs,1) > 1,:);
ioverlap    = ismember(Betas(:,{'UnID','Date'}), overlap);
% Calculate weekly
year   = fix(double(Betas.Date(ikeep & ~ioverlap))/1e4);
week   = weeknum(yyyymmdd2serial(Betas.Date(ikeep & ~ioverlap)))';
[unW,~,subsWeek] = unique([Betas.UnID(ikeep & ~ioverlap) year week],'rows');
dates  = accumarray(subsWeek, Betas.Date(ikeep & ~ioverlap),[size(unW,1),1],@max);
tmp    = accumarray(subsWeek, Betas.Beta(ikeep & ~ioverlap),[size(unW,1),1],@(x) sum(x)/numel(x),NaN);
Betasw = dataset({unW(:,1),'ID'},{dates, 'Date'},{tmp, 'Week'});

clearvars -except Betas Betasd Betasw ikeep resdir
% save debugstate