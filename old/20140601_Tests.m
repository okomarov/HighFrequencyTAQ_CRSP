%% DVD Index file 
d   = 'C:\TAQ\Data\T201101A.IDX';
fid = fopen(d);
data = reshape(fread(fid,inf,'*uint8'),22,[])';
fclose(fid);

out = struct('symbol', cellstr(char(data(:,1:10))),...
    'date'  , num2cell(typecast(reshape(data(:,11:14)',[],1),'uint32')),...
    'start' , num2cell(typecast(reshape(data(:,15:18)',[],1),'uint32')),...
    'end'   , num2cell(typecast(reshape(data(:,19:22)',[],1),'uint32')));

% Transaction file
d   = 'C:\TAQ\Data\T201101A.BIN';
fid = fopen(d);
data = reshape(fread(fid,inf,'*uint8'),19,[])';
fclose(fid);

%% Convert to structure - tested on 1.6GB - 29 sec to save; 10 sec to load; +10% space - [NOT CONVENIENT]
for ii = 1:numel(mst{1})-1
    idx = mst{2}(mst{3}(ii),3):mst{2}(mst{3}(ii+1),3)-1;
    S.(mst{1}{ii}) = {tmp{1}(idx,:) tmp{2}(idx) tmp{3}(idx) tmp{4}(idx,:) tmp{5}(idx,:) tmp{6}(idx)};
end
S.master = mst;
save('C:\TAQ\s.mat','-struct','S','-v7.3')

%% Check that imported the right data - [INCORPORATED INTO IMPORT]
d   = 'C:\TAQ\200912_03.csv';
fid = fopen(d);
opt = {'Delimiter',',','HeaderLines',1};

% Skip all till start-2
sten = [45147978 45164572];
N = 1e5;
tic
for ii = 1:floor(sten(1)/N)
    textscan(fid,'%c%*[^\n]',N,opt{:});
end
toc
textscan(fid,'%c%*[^\n]',rem(sten(1),N)-2,opt{:});
% Import till end
cmp = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',diff(sten)+2,opt{:});
fclose(fid);
% Compare
isequal(cat(2,cmp{3:5}), data{1}(sten(1)-1:sten(2),:))
isequal(         cmp{6}, data{2}(sten(1)-1:sten(2),:))
isequal(         cmp{7}, data{3}(sten(1)-1:sten(2),:))
isequal(cat(2,cmp{8:9}), data{4}(sten(1)-1:sten(2),:))
isequal(  char(cmp{10}), data{5}(sten(1)-1:sten(2),:))
isequal(  char(cmp{11}), data{6}(sten(1)-1:sten(2),:))

[~,~,id] = unique(cmp{1});
blocks = mst{2}(find(sten(1) == mst{2}(:,3))-1:find(sten(2) == mst{2}(:,3)),1:2);
isequal(unique([id cmp{2}],'rows'),[blocks(:,1)-blocks(1)+1 blocks(:,2)])
%% Load speed test

% SEQUENTIAL
% ~3253 sec or 54m
d  = 'E:\data\';
dd = dir(fullfile(d,'*.mat'));

% Load time test (144sec 1:
t = zeros(numel(dd),1);
for ii = 1:1815
    tic
    load(fullfile(d,dd(ii).name))
    t(ii) = toc;
end

% PARALLEL - SLOWER! (usb 2.0 or worker-client data-transfer bottleneck?)
% ~4439 sec or 1h:13m 
% RAM: 340MB matlab + 130MB x worker 
% parallel.cluster.Local('NumWorkers',6) INEFFICIENT
d  = 'E:\data\';
dd = dir(fullfile(d,'*.mat'));
t2 = zeros(numel(dd)-1,1);
%fsize = [dd.bytes]/1048576;
%pos = find(fsize > 90 & fsize < 120);

tic
matlabpool 4
toc
tic
parfor ii = 1:numel(dd)-1
    tic
    s = load(fullfile(d,dd(ii).name));
    s = 1;
    %s = load(fullfile(d,dd(pos(ii)).name));
    t2(ii) = toc;
end
toc
tic
matlabpool close
toc

[srt,idx] = sort([dd(:).bytes]);
plot(srt, t(idx))

%% Save format
d  = 'E:\data\';
tic;load(fullfile(d,'199908_01.mat'));toc
tic;save('test.mat','data','mst','-v6');toc
tic;load('test.mat');toc
%% Textscan test
n   = 50;
opt = {'Delimiter',',','HeaderLines',1};
N   = [500 1e3 5e3 1e4:1e4:1e5 2e5:1e5:1e6];
t   = zeros(n,numel(N));
for jj = 1:23;
    disp(N(jj))
    fid = fopen('C:\TAQ\8e1e9fb052f2b2b6.csv');
    for ii = 1:n
        tic
        foo   = textscan(fid,'%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c',N(jj),opt{:});
        t(ii,jj) = toc;
    end
    fclose(fid);
end
plot(mean(t)*N(end)./N)
set(gca,'Xtick',[2,8 12 22],'XtickL',N([2,8, 12, 22]))
%% Check zips
dd  = 'E:\data';
d   = dir(fullfile(dd,'raw','*.zip'));
idx = false(size(d));
for f = 1:numel(d)
    try
        zipJavaFile  = java.io.File(fullfile(dd,'raw',d(f).name));
        zipFile = org.apache.tools.zip.ZipFile(zipJavaFile);
    catch
        idx(f) = true;
    end
    zipFile.close
end
nnz(idx)
find(idx)
%% Selection 
d = '.\data\TAQ';
s = load(fullfile(d,'T0025.mat'));
% Set to NaN prices that were subject to
% Corrections: all apart 0
% Anomalous sale conditions: all but 'E' (69), 'F' (70) or ' ' (32)
% Zero prices
tmp = uint8(upper(s.data{5}));
pos = s.data{4}(:,2) ~= 0 | ~all( tmp == 32 | tmp == 69 | tmp == 70, 2) |...
      s.data{2} == 0;
s.data{2}(pos) = NaN;

%% Check data.Condition
cd C:\HFbetas 
addpath .\utils
writeto = '.\results\';

% Open matlabpool
if matlabpool('size') == 0, matlabpool 4, end

fun = 'Check_dataCondition';
setupemail
try
    tic
    d   = '.\data\TAQ';
    dd  = dir(fullfile(d,'*.mat'));
    N   = numel(dd);
    res = deal(cell(N,1));
    % LOOP in parallel
    parfor f = 1:N
        disp(f)
        s      = load(fullfile(d,dd(f).name),'data');
        res{f} = unique(s.data.Condition);
    end
     % Collect all results and convert to dataset
    res = cat(1,res{:});
    % Export results and notify
    save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),fun)), 'res')
    message = sprintf('Task ''%s'' terminated in %s',fun,sec2time(toc)); disp(message)
    sendmail('o.komarov11@imperial.ac.uk', message,'')
catch err
    filename = fullfile(writeto, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
    save(filename,'err')
    sendmail('o.komarov11@imperial.ac.uk',sprintf('ERROR in task ''%s''',fun), err.message, {filename})
    rethrow(err)
end
matlabpool close
rmpref('Internet','SMTP_Password')
%% Test return conversion
d = '.\data\TAQ';
s = load(fullfile(d,'T0025.mat'));

% Set to NaN prices that were subject to
% Corrections: all apart 0
% Anomalous sale conditions: all but 'E' (69), 'F' (70) or ' ' (32)
% Zero prices
tmp = uint8(upper(s.data{5}));
pos = s.data{4}(:,2) ~= 0 | ~all( tmp == 32 | tmp == 69 | tmp == 70, 2) |...
      s.data{2} == 0;
s.data{2}(pos) = NaN;

% Set to NaN overnight returns
ret  = s.data{2}(2:end)./s.data{2}(1:end-1)-1;
ret(setdiff(s.mst{2}(:,3)-1,[-1,0])) = NaN;
% Get stats
st = AllStats(ret);

pos = find(ret == st.Min)+1;
scrollto('ret', [pos,1])
scrollto('s.data{2}', [pos,1])
scrollto('s.data{1}', [pos,1])
scrollto('tmp', [pos,1])


pos2 = find(s.mst{2}(:,3) <= pos & s.mst{2}(:,4) >= pos);
scrollto('s.mst{2}',pos2)


% Flag erased rets
idx = false(size(ret));
for r = 1:size(s.mst{2},1)
    fromto     = s.mst{2}(r,3):s.mst{2}(r,4)-1;
    [srt,psrt] = sort(ret(fromto));
end

% Get stats
AllStats(ret)

[~, pmin] = min(ret);
scrollto('s.data{1, 2}',pmin)
scrollto('ret',pmin)
scrollto('tmp',pmin)
pmst = find(s.mst{2}(:,3) < pmin & s.mst{2}(:,4) >= pmin)
scrollto('s.mst{2}',pmst)
scrollto('s.mst{1}',s.mst{2}(pmst,1))

fid = fopen('C:\TAQ\b045446dd7078060.csv');
A = textscan(fid,'%s%f%s%f%*[^\n]','Delimiter',',','HeaderLines',1);
fclose(fid)
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

%% Selection/filtering

% Load big master file
d = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')

% Results directory
resdir = '.\results\';

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
load(fullfile(resdir,'taq2crsp.mat'))

% Loop
% Map dates of the mst to datef
% Map ids to taq2crsp

totbad      = accumarray(mst.Id, mst.Baddays);
totobs      = accumarray(mst.Id,  mst.To - mst.From +1);
% hist(totbad./totobs,100) 
badseries   = totbad./totobs > .1;
mst.Baddays = mst.Baddays | badseries(mst.Id);
clearvars -except mst d SP500

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
perfew       = accumarray(mst.Id, ifewtrades)./accumarray(mst.Id, 1) > .5;
mst.Timestep = ifewtrades | ismember(mst.Id, find(perfew));

%% Betas
cd C:\HFbetas
addpath '.\utils' '.\utils\nth_element'

% Load SPX (tickwrite)
d         = '.\data\';
filename  = unzip(fullfile(d,'Tickwrite','SP.zip'), fullfile(d,'Tickwrite'));
SP500tick = dataset('File',filename{:},'Delimiter',',','ReadVarNames',1,'format','%f%f%f%*[^\n]');
SP500tick = replacedata(SP500tick,@(x) yyyymmdd2serial(x + ((x>9e5).*19e6 + (x<=9e5).*20e6)),'Date');
SP500tick = replacedata(SP500tick,@(x) hhmmssfff2serial(x),'Time');
SP500tick.Datetime = SP500tick.Date + SP500tick.Time;
delete(filename{:})

% Median price for same timestamp
[SP500new,~,subs] = unique(SP500tick(:,{'Date','Datetime'}));
SP500new.Price = accumarray(subs,SP500tick.Price,[],@fast_median);
clear SP500tick subs

% Sample SPX
grid  = (9.5/24:5/(60*24):16/24)';
SP500 = fixedsampling(double([SP500new.Datetime, SP500new.Price]), 'Previous', grid);
SP500 = mat2dataset(SP500,'VarNames',{'Datetime','Price'});

% Load MASTER file
d   = '.\data\TAQ';
load(fullfile(d,'master'),'-mat')

% Prepare pre-cached SP500
tic
mstdates      = unique(mst(:,{'Date','File'}));
mstdates.Date = yyyymmdd2serial(double(mstdates.Date));
toc
N = max(mstdates.File);
tmp = cell(N,1);
for ii = 1:N
    disp(ii)
    dt = mstdates.Date(ii == mstdates.File);
    tmp{ii} = SP500(ismembc(fix(SP500.Datetime),dt), :);
%     tmp{ii}.File = repmat(uint16(ii), size(tmp{ii},1),1);
end
toc
SP500 = tmp;
clear tmp SP500new

% Load daily stats
% Analyze('dailystats',{'Min','Max','MedPrice','Nrets'},[],1)
load .\results\20131230_1732_dailystats.mat
mst = [mst, res];

% Load avg timestes:
% - Worst case 13 trades with an AVERAGE of 30min timestep
% clear res 
% Analyze('avgtimestep','Timestep', mst(:, {'File','MedPrice'}),1)
load .\results\20131230_1844_avgtimestep.mat
ifewtrades   = isnan(res.Timestep) | res.Timestep > 1/48 | mst.Nrets < 12;
perfew       = accumarray(mst.Id, ifewtrades)./accumarray(mst.Id, 1) > .5;
mst.Timestep = ifewtrades | ismember(mst.Id, find(perfew));

% Bad prices days/series
% Analyze('badprices',[],mst(:, {'File','MedPrice'}))
load .\results\20131230_2142_badprices.mat
mst.Baddays = res;
totbad      = accumarray(mst.Id, mst.Baddays);
totobs      = accumarray(mst.Id,  mst.To - mst.From +1);
    % hist(totbad./totobs,100) 
badseries   = totbad./totobs > .1;
mst.Baddays = mst.Baddays | badseries(mst.Id);
clearvars -except mst d SP500

% Run min/max diagnostics
mst = mst(:, {'File','MedPrice','Timestep','Baddays'});
% Analyze('cleanprices',{'Min','Max'},mst,1)

% Cache mst into file blocks
cached = num2cell([accumarray(mst.File,(1:size(mst))',[],@(x) {mst(x,2:end)}), SP500],2);
clear mst SP500
save test

load test
Analyze('betas','Beta',cached)
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