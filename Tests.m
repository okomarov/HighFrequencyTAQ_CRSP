%% Selection/filtering

% Load big master file
d = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')

% Results directory
resdir = '.\results\';

% Map unique ID to mst
testname = 'uniqueID';
try
    dd    = dir(fullfile(resdir, sprintf('*%s.mat', testname)));
    names = sort({dd.name});
    load(fullfile(resdir,names{end}))
catch
    % Re-create mapping
    load(fullfile(resdir,'taq2crsp.mat'))
    taq2crsp = replacedata(taq2crsp,@uint32,'ID');
    % Preallocate
    unID = zeros(size(mst,1),1,'uint16');
    % LOOP for each symbol in mst
    for ii = 1:numel(ids)
        symbol   = ids{ii};
        % Extract records from taq2crsp corresponding to TAQ's symbol
        tmp      = sortrows(taq2crsp(strcmpi(symbol,taq2crsp.symbol), {'ID','datef'}),'datef');
        if isempty(tmp)
            tmp = sortrows(taq2crsp(strcmpi(regexprep(symbol,'p','PR'), taq2crsp.symbol), {'ID','datef'}),'datef');
        end
        if isempty(tmp),fprintf('%d\n',ii),continue,end
        imst     = find(mst.Id == ii);
        % Find to which intervals the records belong
        [~,itmp] = histc(mst.Date(imst), [tmp.datef; 99999999]);
        nnzero   = itmp ~= 0;
        % Assgin unique ID to mst
        unID(imst(nnzero)) = tmp.ID(itmp(nnzero));
    end
    % Unmatched
    mst.UnID(mst.UnID == 0) = intmax('uint16');
    res = dataset(unID,'VarNames','UnID');
    save(fullfile(resdir, sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),testname)), 'res')
end
mst = [mst, res];

% clearvars -except mst d ids resdir
% save debugstate 

% Median and other dailystats
testname = 'dailystats';
try 
    dd    = dir(fullfile(resdir, sprintf('*%s.mat', testname)));
    names = sort({dd.name});
    load(fullfile(resdir,names{end}))
catch
    res = Analyze(testname,{'Min','Max','MedPrice','Nrets'});
end
mst = [mst, res];

% Bad prices days 
testname = 'badprices';
try 
    dd    = dir(fullfile(resdir, sprintf('*%s.mat', testname)));
    names = sort({dd.name});
    load(fullfile(resdir,names{end}))
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
    dd    = dir(fullfile(resdir, sprintf('*%s.mat', testname)));
    names = sort({dd.name});
    load(fullfile(resdir,names{end}))
catch
    res = Analyze(testname,'Timestep', mst(:, {'File','MedPrice'}));
end
mst = [mst, res];

% Select on basis of minimum number of observations
% - Worst case 13 trades with an AVERAGE of 30min timestep
ifewtrades   = isnan(res.Timestep) | res.Timestep > 1/48 | mst.Nrets < 12;
perfew       = accumarray(mst.UnID, ifewtrades)./accumarray(mst.UnID, 1) > .5;
mst.Timestep = ifewtrades | ismember(mst.UnID, find(perfew));

% Sample at 5 min
mst = mst(:, {'File','UnID','MedPrice','Baddays','Timestep'});
testname = 'sample';
Analyze(testname,{'ID','Datetime','Price'}, mst(:, {'File','UnID','MedPrice','Baddays','Timestep'}));

%% Betas
addpath '.\utils' '.\utils\nth_element'
resdir = '.\results';
d      = '.\data\';

% Load SPX (tickwrite)
testname = 'SP500';
try 
    dd    = dir(fullfile(resdir, sprintf('*%s.mat', testname)));
    names = sort({dd.name});
    load(fullfile(resdir,names{end}))
catch
    
    filename  = unzip(fullfile(d,'Tickwrite','SP.zip'), fullfile(d,'Tickwrite'));
    SP500tick = dataset('File',filename{:},'Delimiter',',','ReadVarNames',1,'format','%f%f%f%*[^\n]');
    SP500tick = replacedata(SP500tick,@(x) yyyymmdd2serial(x + ((x>9e5).*19e6 + (x<=9e5).*20e6)),'Date');
    SP500tick = replacedata(SP500tick,@(x) hhmmssfff2serial(x),'Time');
    SP500tick.Datetime = SP500tick.Date + SP500tick.Time;
    delete(filename{:})
    
    % Median price for same timestamp
    [SP500,~,subsID] = unique(SP500tick(:,{'Date','Datetime'}));
    SP500.Price = accumarray(subsID,SP500tick.Price,[],@fast_median);
    save(fullfile(resdir,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),testname)), 'SP500')
end

% Sample SP500 (until 3:15 PM!)
grid  = (9.5/24:5/(60*24):(15+15/60)/24)';
ngrid = numel(grid);
SP500 = fixedsampling(double([SP500.Datetime, SP500.Price]), 'Previous', grid);
SP500 = mat2dataset(SP500,'VarNames',{'Datetime','Price'});

% SP500 ret (zeroing overnight)
SPret = [SP500.Datetime(2:end) SP500.Price(2:end)./SP500.Price(1:end-1)-1];
SPret = SPret(diff(fix(SP500.Datetime)) == 0,:);

% Cache SP500 returns by days
SPdays = unique(fix(SPret(:,1)),'stable');
SPret  = mat2cell(SPret(:,2),repmat(ngrid-1,size(SPret,1)/(ngrid-1),1));


% LOOP by sampled data 
d = '.\data\TAQ\sampled';
dd = dir(fullfile(d,'S*.mat'));
res = cell(numel(dd),1);
parpool(4)
parfor f = 1:numel(dd)
    disp(f)
    % Load data and get returns until 3:15 PM
    s      = load(fullfile(d,dd(f).name));
    ret    = [s.data.Datetime(2:end), s.data.Price(2:end)./s.data.Price(1:end-1)-1];
    idx    = rem(ret(:,1),1) <= (15+16/60)/24;
    ret    = ret(idx,:);
    ret    = ret([true; diff(rem(ret(:,1),1)) >= 0],:);
    
    % Map SP500 rets to stock rets
    days    = yyyymmdd2serial(double(s.mst.Date));
    spret   = cat(1,SPret{ismembc2(days, SPdays)});
    prodret = spret.*ret(:,2);
    subsID    = reshape(repmat(1:size(s.mst,1),ngrid-1,1),[],1);
    ikeep   = ~isnan(prodret);
    beta    = accumarray(subsID(ikeep), prodret(ikeep))./accumarray(subsID(ikeep), spret(ikeep).^2);
    
    % Store results
    res{f} = s.mst(:,{'UnID','Date'});
    res{f}.Beta = beta;
end
Betas = cat(1,res{:});
save(fullfile(resdir,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'Betas')), 'Betas')
%% Smooth Betas
resdir = '.\results';

% Load Betas
dd     = dir(fullfile(resdir, sprintf('*%s.mat', 'Betas')));
names  = sort({dd.name});
load(fullfile(resdir,names{end}))

% Sort betas
Betas = sortrows(Betas,{'UnID','Date'});

% Index by ID
[unID, ~,subsID] = unique(Betas.UnID);

% Exclude series with less than 20 observations, i.e. one month of data
ifewdays = accumarray(subsID,1) < 20;
ikeep    = ~ismember(Betas.UnID, unID(ifewdays));

% Moving averages
% tmp = accumarray(subsID(ikeep), Betas.Date(ikeep),[],@issorted, true);
sz      = size(unID);
Res.SMA = accumarray(subsID(ikeep), Betas.Beta(ikeep),sz,@(x) {conv(x,ones(1,5)/5,'valid')});
Res.EMA = accumarray(subsID(ikeep), Betas.Beta(ikeep),sz,@(x) {movavg(x,1,5,'e')});

% Weekly betas
[unW,~,subsWeek] = unique([Betas.UnID fix(double(Betas.Date)/1e4) weeknum(yyyymmdd2serial(double(Betas.Date)))'],'rows');
Res.Week         = accumarray(subsWeek(ikeep), Betas.Beta(ikeep),[size(unW,1),1],@(x) sum(x)/numel(x),NaN);

%% Size quantiles
resdir = '.\results';
vars = {'cusip','symbol','datef'};
load(fullfile(resdir, 'taq2crsp.mat'))

% Load shrout getting rid of 0s
load(fullfile(resdir, 'TAQshrout.mat'))
TAQshrout.Properties.VarNames = [vars, 'shrout'];

% Direct match
TAQshrout.ID      = zeros(size(TAQshrout,1),1,'uint16');
[idx,pos]         = ismember(TAQshrout(:,vars(1:2)), taq2crsp(:,vars(1:2)));
TAQshrout.ID(idx) = taq2crsp.ID(pos(idx));

% Organize number of shares as monthly panel





%% Verify MANUAL vs AUTOMATIC betas
addpath .\utils\ .\utils\nth_element\ .\utils\MFE

% Load automatic betas
load .\data\TAQ\master -mat
load .\results\20131231_1338_betas.mat
mst.Date = yyyymmdd2serial(double(mst.Date));
idx      = mst.Date >= 730124 &...
           mst.Date <= 730485 &...
           mst.Id == find(strcmpi(ids,'AAPL'));
autobetas = res.Beta(idx);
clear res mst

% Load AAPL sample from 1999
filename = unzip('.\data\TAQ\datachk\1999\AAPL_csv.zip');
fid        = fopen(filename{end});
data       = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c','Delimiter',',','HeaderLines',1);
data{1,10} = char(data{1,10});
szCond     = size(data{1,10});
data{1,10} = [data{1,10} repmat(' ',szCond(1),2-szCond(2))];
data       = [data(2) cat(2,data{3:5}) data(6:7) cat(2,data{8:9}) data(10:11)];
fclose(fid);
delete(filename{1});
% Into dataset
data = dataset({yyyymmdd2serial(double(data{1})), 'Date'},...
                { yyyymmdd2serial(double(data{1})) ...
                + hhmmssmat2serial(data{2}), 'Datetime'},...
                {data{3},'Price'},...
                {data{4},'Volume'},...
                {data{5},'G127_Correction'},...
                {data{6},'Condition'});
% Select trades
data = data(~selecttrades(data),:);
load .\results\20131230_1732_dailystats.mat
mst  = res(idx,:);
nobs = diff([0; find(diff(data.Date)); size(data,1)]);
[~,ibadprice] = histc(data.Price./RunLength(mst.MedPrice, nobs), [.65,1.51]);
ibadprice     = ibadprice ~= 1;
data(ibadprice,:) = [];
     
% Load SP 
d         = '.\data\';
filename  = unzip(fullfile(d,'Tickwrite','SP.zip'), fullfile(d,'Tickwrite'));
SP500tick = dataset('File',filename{:},'Delimiter',',','ReadVarNames',1,'format','%f%f%f%*[^\n]');
SP500tick = replacedata(SP500tick,@(x) yyyymmdd2serial(x + ((x>9e5).*19e6 + (x<=9e5).*20e6)),'Date');
SP500tick = replacedata(SP500tick,@(x) hhmmssfff2serial(x),'Time');
SP500tick.Datetime = SP500tick.Date + SP500tick.Time;
delete(filename{:})

idx = SP500tick.Date >= fix(data.Datetime(1)) & SP500tick.Date <= fix(data.Datetime(end));
SP500tick = SP500tick(idx,:);

save .\results\checkBetasData SP500tick data autobetas

load .\results\checkBetasData

% Median price for same timestamp
[SP500new,~,subs] = unique(SP500tick(:,{'Date','Datetime'}));
SP500new.Price    = accumarray(subs,SP500tick.Price,[],@fast_median);

[datanew,~,subs] = unique(data(:,{'Date','Datetime'}));
datanew.Price    = accumarray(subs,data.Price,[],@fast_median);

% Sample
grid  = (9.5/24:5/(60*24):16/24)';
SP500sampled = fixedsampling([SP500new.Datetime, double(SP500new.Price)], 'Previous', grid);
SP500sampled = mat2dataset(SP500sampled,'VarNames',{'Datetime','Price'});

datasampled = fixedsampling([datanew.Datetime, double(datanew.Price)], 'Previous', grid);
datasampled = mat2dataset(datasampled,'VarNames',{'Datetime','Price'});

% Returns
len                   = length(grid);
dataret               = datasampled.Price(2:end)./datasampled.Price(1:end-1)-1;
dataret(len:len:end)  = [];
SP500ret              = SP500sampled.Price(2:end)./SP500sampled.Price(1:end-1)-1;
SP500ret(len:len:end) = [];

% Beta
subs       = repmat(1:length(dataret)/(len-1),len-1,1);
rcmanual   = accumarray(subs(:),dataret.*SP500ret,[],@nansum);
vmanual    = accumarray(subs(:),SP500ret.^2,[],@nansum);
betamanual = rcmanual./vmanual;
             

% Compare against MFE toolbox by Sheppard
undays = unique(SP500new.Date)';
ndays  = numel(undays);
[rc,v] = deal(zeros(ndays,1));
type = 'fixed';%'CalendarTime';
ggrid = grid;%5/(60*24);

nsamples = 1;
for ii = 1:ndays
    idata  = undays(ii) == datanew.Date;
    isp500 = undays(ii) == SP500new.Date;
    if ~isscalar(ggrid)
        ggrid = undays(ii) + grid;
    end
    
    tmp    = realized_covariance( datanew.Price(idata), datanew.Datetime(idata), ...
                                 SP500new.Price(isp500), SP500new.Datetime(isp500),...
                                 'unit',type,ggrid,nsamples);
    rc(ii) = tmp(2);
    v(ii)  = realized_variance(SP500new.Price(isp500), SP500new.Datetime(isp500),...
                                 'unit',type,ggrid,nsamples);
end
betaSheppard = rc./v;

plot([betamanual,  autobetas,betaSheppard]), legend('manual','auto','mfe')
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

