%% Counts selection rule
path2data = '.\data\TAQ';
testname  = 'selrulecounts';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,[],[],fullfile(path2data,'T*.mat'));
end

% Plot 
for ii = 1:size(res,2)-1
    if ii == 3
        data          = cat(1,res{:,ii});
        idx           = all(data.Val == ' ' | data.Val == 'E' | data.Val == 'F' | data.Val == '@',2);
        data.Val      = cellstr(data.Val); 
        data.Val(idx) = {'  '};
        data = unstack(data,'Count','Val');
    else
        data = unstack(cat(1,res{:,ii}),'Count','Val');
    end
    vnames = getVariableNames(data);
    dates  = datenum(double(data.yyyymm/100), double(rem(data.yyyymm,100)+1), 1)-1;
    data   = table2array(data(:,2:end));
    data(isnan(data)) = 0;
    
    figure
    subplot(211)
    title('Montly counts of price observations by type (absolute and %)')
    area(dates, data)
    dynamicDateTicks, axis tight
    
    subplot(212)
    area(dates, bsxfun(@rdivide, data, sum(data,2))*100)
    dynamicDateTicks, axis tight
    legend(vnames(2:end),'Location','west');
end
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
    loadresults('TAQmaster')
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
loadresults('TAQmaster')
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
loadresults('taq2crsp')
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
    loadresults('msenames')
catch
    msenames = importMsenames('.\data\CRSP\');
end
msenames = unique(msenames(:,{'PERMNO','NAMEDT','SHRCD'}));
% Time consolidation
idx      = isfeatchange(msenames(:,[1,3,2]));                 
msenames = msenames(idx,:);

% Load taq2crsp
loadresults('taq2crsp')
taq2crsp       = taq2crsp(~isnan(taq2crsp.permno),:);

% Sort (match by permno)
taq2crsp       = sortrows(taq2crsp,{'permno','datef'});
taq2crsp.Shrcd = zeros(size(taq2crsp,1),1,'uint8');

% Add SHRCD to TAQ2CRSP
% -------------------------------------------------------------------------
% Single entry link between msenames and TAQ2CRSP
idx = [true; logical(diff(msenames.PERMNO))];
for ii = 2:numel(idx)
    if ~idx(ii)
        idx(ii-1) = false;
    end
end
[~,posAll] = ismember(taq2crsp.permno, msenames.PERMNO);
idirect    = ismember(posAll, find(idx));
taq2crsp.Shrcd(idirect) = msenames.SHRCD(posAll(idirect));

% Multiple entry link
for ii = 1:numel(idx)
    if ~idx(ii)
        permno = msenames.PERMNO(ii);
        to     = find(msenames.PERMNO == permno,1, 'last');
        dates  = [msenames.NAMEDT(ii:to); inf];
        shrcd  = [0; msenames.SHRCD(ii:to)];
        
        itaq                 = taq2crsp.permno == permno;
        [~,datebin]          = histc(taq2crsp.datef(itaq), dates);
        taq2crsp.Shrcd(itaq) = shrcd(datebin+1);
    end
end

% Add SHRCD to MASTER
% -------------------------------------------------------------------------
% Re-sort (match by symbol)
taq2crsp = unique(taq2crsp(:,{'symbol','datef','Shrcd'}));
% Time consolidation
idx      = isfeatchange(taq2crsp(:,{'symbol','Shrcd','datef'}),2:3);
taq2crsp = taq2crsp(idx,:);

% Master file
path2data  = '.\data\TAQ';
load(fullfile(path2data, 'master'), '-mat');

% Preallocation
ids = regexprep(ids,'p','PR');
ids = regexprep(ids,'\.','');
Shrcd = zeros(size(mst,1),1,'uint8');

% p = poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'}, 'debug',false);

% Link to CRSP ~1hr
[unsymbols,~,subsymbols] = unique(taq2crsp.symbol);
tic
for ii = 1:numel(unsymbols)
    disp(ii)
    % taq2crsp info
    symbol = unsymbols{ii};

    % master info
    id = find(strcmpi(ids,symbol));
    if isempty(id)
        fprintf('No match: %s - iter %d\n', symbol,ii), continue
    end
    if numel(id) == 1
        imst = mst.Id == id;
    else
        imst = ismember(mst.Id,id);
    end
    
    itaq   = subsymbols == ii;
    dates  = [taq2crsp.datef(itaq);inf];
    shrcds = uint32([0; taq2crsp.Shrcd(itaq)]);
    
    [~, datebin] = histc(mst.Date(imst), dates);
    Shrcd(imst)  = shrcds(datebin+1);
end
toc

res = mst(:,{'Id','Date'});
res.Shrcd = Shrcd;

% Save
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'shrcd')), 'res')

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
counts   = dataset2table(unstack(counts,'Counts','Shrcd'));
counts.Properties.VariableNames{1} = 'Dates';
vnames = getVariableNames(counts);

% Save
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'shrcdcounts')), 'counts')

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
%% Selection/filtering

% Load big master file
d = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')

% Map unique ID to mst
testname = 'uniqueID';
try
    loadresults(testname,'res')
% Re-create mapping
catch
    tic
    res = mapUnid2mst(mst, ids);
    sec2time(toc)
end
[~,pos] = ismember(mst(:,{'Id','Date'}), res(:,{'Id','Date'}));
% clearvars -except mst d ids resdir

% Median and other dailystats
testname = 'dailystats';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,{'Min','Max','MedPrice','Nrets'});
end
mst = [mst, res];

% Bad prices days
testname = 'badprices';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,'Baddays',mst(:, {'File','MedPrice'}));
end
mst = [mst, res];

% Bad series
totbad      = accumarray(mst.UnID, mst.Baddays);
totobs      = accumarray(mst.UnID, mst.To - mst.From +1);
% hist(totbad./totobs,100)
badseries   = totbad./totobs > .1;
badseries(end) = true; % for the unmatched
mst.Baddays = mst.Baddays | badseries(mst.UnID);

% clearvars -except mst d ids resdir
% mst(:,{'Min','Max'}) = [];

% Average time step
testname = 'avgtimestep';
try
    loadresults(testname,'res')
catch
    res = Analyze(testname,'Timestep', mst(:, {'File','MedPrice'}));
end
mst = [mst, res];

% Select on basis of minimum number of observations
% - Worst case 13 trades with an AVERAGE of 30min timestep
ifewtrades   = isnan(res.Timestep) | res.Timestep > 1/48 | mst.Nrets < 12;
perfew       = accumarray(mst.UnID, ifewtrades)./accumarray(mst.UnID, 1) > .5;
mst.Timestep = ifewtrades | perfew(mst.UnID);

% Sample at 5 min
mst = mst(:, {'File','Id','UnID','MedPrice','Baddays','Timestep'});
testname = 'sample';
Analyze(testname,{'ID','Datetime','Price'}, mst);
%% Count 0 returns

testname = 'countnullrets';
try
    loadresults(testname,'counts')
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
    loadresults(testname,'counts')
    if isa(counts,'dataset'), counts = dataset2table(counts); end
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
mktcap = getMktCap(unids, refdates);

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
%% Betas from SP500proxy
addpath '.\utils' '.\utils\nth_element'
testname = 'betas';
try
    loadresults(testname)
catch
    path2data = '.\data\TAQ\sampled';
    
    % Load SPY (etf)
    loadname = 'sp500proxy';
    try
        loadresults(loadname, 'spproxy')
    catch
        spproxy = sp500intraday;
    end
    
    % SP500 ret (zeroing overnight)
    spret = [spproxy.Datetime(2:end) spproxy.Price(2:end)./spproxy.Price(1:end-1)-1];
    spret = spret(diff(rem(spproxy.Datetime,1)) >= 0,:);
    
    % Cache SP500 returns by days
    load(fullfile(path2data,'master'),'-mat','mst');
    nfiles   = max(mst.File);
    cached   = cell(nfiles,1);
    
    dates  = fix(spret(:,1));
    [spdays,~,subs] = unique(dates,'stable');
    spret = mat2cell(spret(:,2), accumarray(subs,1),1);
    
    unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(unique(x))});
    for ii = 1:nfiles
        pos        = ismembc2(unMstDates{ii}, spdays);
        nnzero     = pos ~= 0;
        isp       = ismembc(spdays, unMstDates{ii});
        cached{ii} = {spret(isp) spdays(pos(nnzero))};
    end
    
    % Calculate betas
    Betas = Analyze(testname, [], cached, fullfile(path2data,'S5m_*.mat'));
end
%% Betas from SPY
addpath '.\utils' '.\utils\nth_element'
resdir    = '.\results';

testname = 'betas';
try
    loadresults(testname)
catch
    path2data = '.\data\TAQ\sampled';
    
    % Load SPY (etf)
    loadname = 'spysampled';
    try
        loadresults(loadname, 'SPY')
    catch
        master    = load(fullfile(path2data,'master'),'-mat');
        SPY       = getTaqData(master, 'SPY',[],[],'Price',path2data);
        save(fullfile(resdir,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),loadname)), 'spy')
    end
    
    % SP500 ret (zeroing overnight)
    spyret = [SPY.Datetime(2:end) SPY.Price(2:end)./SPY.Price(1:end-1)-1];
    spyret = spyret(diff(rem(SPY.Datetime,1)) >= 0,:);
    
    % Cache SP500 returns by days
    load(fullfile(path2data,'master'),'-mat','mst');
    nfiles   = max(mst.File);
    cached   = cell(nfiles,1);
    
    dates  = fix(spyret(:,1));
    [spdays,~,subs] = unique(dates,'stable');
    spyret = mat2cell(spyret(:,2), accumarray(subs,1),1);
    
    unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(unique(x))});
    for ii = 1:nfiles
        pos        = ismembc2(unMstDates{ii}, spdays);
        nnzero     = pos ~= 0;
        ispy       = ismembc(spdays, unMstDates{ii});
        cached{ii} = {spyret(ispy) spdays(pos(nnzero))};
    end
    
    % Calculate betas
    Betas = Analyze(testname, [], cached, fullfile(path2data,'S5m_*.mat'));
end

%% Beta quantiles
addpath .\utils\

% Filter settings
sp500only  = true;
commononly = true;

% Get Betas
Betas = getBetas(sp500only, commononly);

% Sample/expand
refdates = serial2yyyymmdd(datenum(1993,2:234,1)-1);
Betas    = sampledates(Betas,refdates,1);

% All days
% refdates = Betas.Date;
% tmp      = table2array(Betas);

% Plot
plotdates = datetime(yyyymmdd2serial(refdates),'ConvertFrom','datenum');
plot(plotdates, prctile(Betas(:,2:end),10:10:90,2))
legend(arrayfun(@(x) sprintf('%d^{th} ',x),10:10:90,'un',0))

if sp500only
    title 'Cross-sectional percentiles of SP500 Betas (with proxy)'
    saveas(gcf, '.\results\BetasSP500_proxy.png')
else
    title 'Cross-sectional percentiles of all Betas (with proxy)'
    saveas(gcf, '.\results\BetasAll_proxy.png')
end
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
mktcap = getMktCap(unids, refdates);

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
mktcap = getMktCap(unids, refdates);

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
    loadresults(loadname, 'spy')
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
loadresults('sp500proxy', 'spproxy')

% Betas*proxy
spproxy      = iprice2dret(spproxy);
[~,pos]      = ismember(Betas.Date,spproxy.Date);
ret          = [NaN; spproxy.Ret];
Betas.Sysret = ret(pos+1).*Betas.Beta;

% Daily rets
try
    loadresults('dailyret','rets')
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
loadresults('spysampled','SPY')

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
loadresults('taq2crsp')
loadresults('betas','betas2')
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
loadresults('betas','Betas')

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