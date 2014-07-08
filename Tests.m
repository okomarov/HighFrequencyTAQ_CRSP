%% Selection rule counts
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
    vnames = data.Properties.VariableNames(2:end);
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
    legend(vnames,'Location','west');
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
sample = getData(master, symbol, date,date);
i40    = sample.G127_Correction(:,1) == 40;
igood  = sample.G127_Correction(:,1) == 0;
x      = hhmmssmat2serial(sample.Time);
plot(x(i40), sample.Price(i40),'xr', x(igood), sample.Price(igood),'.b',...
     x, sample.Price, '-g')

%% Type counts
try
    loadresults('TAQmaster')
catch
    TAQmaster = importMst('.\data\TAQ\raw');
end

TAQmaster = unique(TAQmaster(:, {'SYMBOL','FDATE','TYPE'}));

% Time consolidation
idx        =  [true; ~(strcmpi(TAQmaster.SYMBOL(2:end), TAQmaster.SYMBOL(1:end-1)) & ...
    TAQmaster.TYPE  (2:end)==TAQmaster.TYPE(1:end-1))];
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
typecounts   = dataset2table(unstack(typecounts,'Counts','Val'));
vnames = {'Common', 'Preferred', 'Warrant', 'Right', 'Other', 'Derivative'};
typecounts.Properties.VariableNames = ['Dates' vnames];

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
%% Selection/filtering

% Load big master file
d = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')

% Results directory
resdir = '.\results\';

% Map unique ID to mst
testname = 'uniqueID';
try
    loadresults(testname,'res')
    
    % Re-create mapping
catch
    tic
    % Get taq2crsp ready
    loadresults('taq2crsp')
    taq2crsp = sortrows(taq2crsp,{'symbol','datef'});
    
    % Preallocate
    unID = zeros(size(mst,1),1,'uint16');
    % LOOP for each symbol in mst
    for ii = 1:numel(ids)
        symbol   = ids{ii};
        % Extract records from taq2crsp corresponding to TAQ's symbol
        isymbol  = strcmpi(symbol,taq2crsp.symbol);
        tmp      = taq2crsp(isymbol, {'ID','datef'});
        %         if isempty(tmp)
        %             % Preferred stocks symbol use sometimes the lowercase suffix
        %             % 'p' instead of 'PR' (see daily TAQ guide)
        %             isymbol = strcmpi(regexprep(symbol,'p','PR'), taq2crsp.symbol);
        %             tmp     = taq2crsp(isymbol, {'ID','datef'});
        %         end
        if isempty(tmp),fprintf('%d\n',ii),continue,end
        imst     = find(mst.Id == ii);
        % Find to which intervals the records belong
        [~,itmp] = histc(mst.Date(imst), [tmp.datef; 99999999]);
        nnzero   = itmp ~= 0;
        % Assgin unique ID to mst
        unID(imst(nnzero)) = tmp.ID(itmp(nnzero));
    end
    % Unmatched
    unID(unID == 0) = intmax('uint16');
    res = dataset(unID,'VarNames','UnID');
    save(fullfile(resdir, sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),testname)), 'res')
end
mst = [mst, res];
sec2time(toc)
% clearvars -except mst d ids resdir
% save debugstate

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
%% Betas
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
        SPY       = getData(master, 'SPY',[],[],'Price',path2data);
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
    
    unMstDates = accumarray(mst.File, mst.Date,[],@(x){yyyymmdd2serial(double(unique(x)))});
    for ii = 1:nfiles
        pos        = ismembc2(unMstDates{ii}, spdays);
        nnzero     = pos ~= 0;
        ispy       = ismembc(spdays, unMstDates{ii});
        cached{ii} = {spyret(ispy) spdays(pos(nnzero))};
    end
    
    % Calculate betas
    Betas = Analyze(testname, [], cached, fullfile(path2data,'S5m_*.mat'));
end
%% Smooth Betas
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
week   = weeknum(yyyymmdd2serial(double(Betas.Date(ikeep & ~ioverlap))))';
[unW,~,subsWeek] = unique([Betas.UnID(ikeep & ~ioverlap) year week],'rows');
dates  = accumarray(subsWeek, Betas.Date(ikeep & ~ioverlap),[size(unW,1),1],@max);
tmp    = accumarray(subsWeek, Betas.Beta(ikeep & ~ioverlap),[size(unW,1),1],@(x) sum(x)/numel(x),NaN);
Betasw = dataset({unW(:,1),'ID'},{dates, 'Date'},{tmp, 'Week'});

clearvars -except Betas Betasd Betasw ikeep resdir
% save debugstate
%% Beta quantiles
load debugstate
% DAILY
% Remove overlapping
[un,~,subs] = unique(Betasd(:,{'ID','Date'}));
overlap     = un(accumarray(subs,1) > 1,:);
[~,idx]     = setdiff(Betasd(:,1:2), overlap);
% taq2crsp(ismember(taq2crsp.permno, taq2crsp.permno(taq2crsp.ID == overlap.ID(n))) | taq2crsp.ID == overlap.ID(n),:)

% Pivot
tmp = Pivot([double(Betasd.ID(idx)), double(Betasd.Date(idx)) double(Betasd.SMA(idx))]);

% Add reference dates
refdates       = serial2yyyymmdd(datenum(1993,2:234,1)-1);
alldates       = union(tmp(2:end,1), refdates);
[~,pdates]     = ismember(tmp(2:end,1),alldates);
data           = NaN(numel(alldates),size(tmp,2)-1);
data(pdates,:) = tmp(2:end,2:end);
tmp            = [[NaN; alldates], [tmp(1,2:end); data]];

% Fill in previous value
tmp(2:end,2:end) = nanfillts(tmp(2:end,2:end),1);

% Keep reference dates only
tmp = tmp([true; ismember(tmp(2:end,1), refdates)],:);

% Plot
plot(yyyymmdd2serial(refdates), prctile(tmp(2:end,2:end),10:10:90,2))
dynamicDateTicks
title 'Cross-sectional percentiles of 5-day smoothed (SMA) betas'
legend(arrayfun(@(x) sprintf('%d^{th} ',x),10:10:90,'un',0))

%% Size quantiles
addpath .\utils\

vars = {'cusip','symbol','datef'};
loadresults('taq2crsp')

% Load shrout getting rid of 0s
loadresults('TAQshrout')
TAQshrout.Properties.VarNames = [vars, 'shrout'];

% Filter out matches (max score 42, so 100 takes all)
iscore = taq2crsp.score < 100 | isnan(taq2crsp.score);

% Direct match
TAQshrout.ID      = zeros(size(TAQshrout,1),1,'uint16');
[idx,pos]         = ismember(TAQshrout(:,vars(1:2)), taq2crsp(iscore,vars(1:2)));
IDs               = taq2crsp.ID(iscore);
TAQshrout.ID(idx) = IDs(pos(idx));

% Matched
imatched = TAQshrout.ID ~= 0;

% Exclude overlapping records
[un,~,subs] = unique(TAQshrout(imatched,{'ID','datef'}));
overlap     = un(accumarray(subs,1) > 1,:);
imatched    = imatched & ~ismember(TAQshrout(:, {'ID','datef'}), overlap);
% taq2crsp(ismember(taq2crsp.permno, taq2crsp.permno(taq2crsp.ID == overlap.ID(3))),:)

% Get panel of shrout
shrout = Pivot([double(TAQshrout.ID(imatched)), double(TAQshrout.datef(imatched)), TAQshrout.shrout(imatched)]);

% Monthly dates of interest
refdates = serial2yyyymmdd(datenum(1993,2:234,1)-1);
alldates = union(shrout(2:end,1), refdates);

[~,pdates] = ismember(shrout(2:end,1),alldates);
data       = NaN(numel(alldates),size(shrout,2)-1);
data(pdates,:) = shrout(2:end,2:end);
shrout = [[NaN; alldates], [shrout(1,2:end); data]];

% Fill in previous value
shrout(2:end,2:end) = nanfillts(shrout(2:end,2:end));

% Get reference dates
shrout = shrout([true; ismember(shrout(2:end,1), refdates)],:);
%% SP500 betas
addpath .\utils\
spdir  = '.\data\SP500';

% Load TAQ2CRSP
loadresults('taq2crsp')

% Load SP consituents > 31/12/1992
SPconst = dataset('File',fullfile(spdir,'dsp500list.csv'),'Delimiter',',','ReadVarNames',1);
SPconst = SPconst(SPconst.ending > 19921231,:);

% Load Betas
loadresults('Betas')

% Filter out betas (lowest UnID for each PERMNO)
selected = unique(taq2crsp(ismember(taq2crsp.permno, SPconst.PERMNO),{'permno','ID'}), 'permno','first');
Betas    = Betas(ismember(Betas.UnID, selected.ID),{'UnID','Date','Beta'});
[~,pos]  = ismember(Betas.UnID,selected.ID);
Betas.PERMNO = selected.permno(pos);

% Overlapping Betas
[un,~,subs] = unique(Betas(:,{'PERMNO','Date'}));
overlap     = un(accumarray(subs,1) > 1,:);
Betas       = Betas(~ismember(Betas(:, {'PERMNO','Date'}), overlap),:);
% taq2crsp(ismember(taq2crsp.permno, taq2crsp.permno(taq2crsp.ID == overlap.UnID(1))),:)

% Pivot betas
tmp            = Pivot([double(Betas.PERMNO), double(Betas.Date), Betas.Beta]);
refdates       = serial2yyyymmdd(datenum(1993,2:234,1)-1);
alldates       = union(tmp(2:end,1), refdates);
[~,pdates]     = ismember(tmp(2:end,1),alldates);
data           = NaN(numel(alldates),size(tmp,2)-1);
data(pdates,:) = tmp(2:end,2:end);
tmp            = [[NaN; alldates], [tmp(1,2:end); data]];
% Fill in previous value and get only reference dates
tmp(2:end,2:end) = nanfillts(tmp(2:end,2:end),1);
Betas = tmp([true; ismember(tmp(2:end,1), refdates)],:);

% Pivot SP membership
SPconst        = [[double(SPconst.PERMNO)  double(SPconst.start)   ones(size(SPconst,1),1)];
    [double(SPconst.PERMNO)  double(SPconst.ending) -ones(size(SPconst,1),1)]];
tmp            = Pivot(SPconst,[],[],0);
refdates       = serial2yyyymmdd(datenum(1993,2:234,1)-1);
alldates       = union(tmp(2:end,1), refdates);
[~,pdates]     = ismember(tmp(2:end,1),alldates);
data           = zeros(numel(alldates),size(tmp,2)-1);
data(pdates,:) = tmp(2:end,2:end);
tmp            = [[NaN; alldates], [tmp(1,2:end); data]];
% Fill in previous value and get only reference dates
tmp(2:end,2:end) = cumsum(tmp(2:end,2:end));
SPconst = tmp([true; ismember(tmp(2:end,1), refdates)],:);

% Intersect matrices
[~,ia,ib] = intersect(Betas(1,:),SPconst(1,:));

% NaN out when non members
mask = ~logical(SPconst(2:end, ib));
tmp = Betas(2:end, ia);
tmp(mask) = NaN;

% Plot
plot(yyyymmdd2serial(refdates), prctile(tmp,10:10:90,2))
dynamicDateTicks
title 'Cross-sectional percentiles of un-smoothed SP500 Betas'
legend(arrayfun(@(x) sprintf('%d^{th} ',x),10:10:90,'un',0))
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

if isempty(gcp('nocreate'))
    parpool(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
end

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
        day           = yyyymmdd2serial(double(s.mst.Date(ii)));
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
spys      = getData(struct(master), 'spy', date, date,[], path2data);
aapls     = getData(struct(master), 'aapl', date, date,[], path2data);

% Full
path2data = '.\data\TAQ\';
master    = load(fullfile(path2data, 'master'), '-mat');
spy       = getData(struct(master), 'spy', date, date,[], path2data);
aapl      = getData(struct(master), 'aapl', date, date,[], path2data);

% Select
spy  = spy(~selecttrades(spy),{'Datetime','Price'});
aapl = aapl(~selecttrades(aapl),{'Datetime','Price'});

% Median
[dates, ~,subs] = unique(spy.Datetime);
prices          = accumarray(subs,spy.Price,[],@fast_median);
spy             = table(dates,prices,'VariableNames', spy.Properties.VariableNames);

[dates, ~,subs] = unique(aapl.Datetime);
prices          = accumarray(subs,aapl.Price,[],@fast_median);
aapl            = table(dates,prices,'VariableNames', aapl.Properties.VariableNames);

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

