%% Options
OPT_VW = true;

OPT_LAG    = 1;
OPT_PTF_UN = 10;

%% Data
master = loadresults('master');

myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

% Returns
dsf    = loadresults('dsf');
ret    = myunstack(dsf,'Ret');
date   = ret.Date;
permno = ret.Properties.VariableNames(2:end);
ret    = double(ret{:,2:end});

% End-of-Month ismicro
dsf.IsMicro = isMicrocap(dsf,'Prc');
isMicro     = myunstack(dsf,'IsMicro');
[~,pos]     = unique(isMicro.Date/100,'last');
isMicro     = isMicro{pos,2:end};

% Mkt cap
dsf     = getMktCap(dsf);
cap     = myunstack(dsf,'Cap');
[~,pos] = unique(cap.Date/100,'last');
cap     = cap{pos,2:end};

% Realized skewness
rskew = loadresults('skew');
num   = myunstack(rskew,'Num');
den   = myunstack(rskew,'Den');
rskew = num{:,2:end}./den{:,2:end};
clear den num

% Overnight return
reton = loadresults('reton');

% Beta components
beta              = loadresults('beta5minon');
num               = myunstack(beta,'Num');
den               = myunstack(beta,'Den');
beta              = cat(3,num{:,2:end},den{:,2:end});
beta(isinf(beta)) = NaN; % permno 46288 on 19931008 is delisted with close-to-close return of -100%
clear den num

% Factors
ff = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
ff = ff(ismember(ff.Date, unique(dsf.Date)),:);
%% Signals

% Low freqeuncy signals
[signals_LF, hpr, rf] = make_signals(ret,date,ff);
nsig                  = size(signals_LF,3);
snames                = {'alpha','skewh','skewr','betas'};
correlations          = corrxs(signals_LF,snames);

% High freqeuncy signals
signals_HF = make_signals_HF(xstr2num(permno),date,master,reton,ff,rskew,beta);

inan                             = any(isnan(signals_LF),3) | any(isnan(signals_HF),3);
signals_HF(repmat(inan,[1,1,3])) = NaN;


corrxs(cat(3,signals_LF(:,:,4),signals_HF(:,:,3)))


%% Lag
% End-of-Month
signals_LF = signals_LF(1:end-OPT_LAG,:,:);
signals_HF = signals_HF(1:end-OPT_LAG,:,:);
isMicro    = isMicro(1:end-OPT_LAG,:);
cap        = cap(1:end-OPT_LAG,:);

% Lag forward
hpr = hpr(1+OPT_LAG:end,:);
rf  = rf(1+OPT_LAG:end,:);
%% Filter micro
hpr(isMicro) = NaN;
%% PTFRET
if OPT_VW
    opts = struct('PortfolioNumber',OPT_PTF_UN, 'Weights',cap);
else
    opts = struct('PortfolioNumber',OPT_PTF_UN);
end

% Alpha
[ptfret{1},~,counts{1},avgsig{1}] = portfolio_sort(hpr, signals_LF(:,:,1), opts);

% Alpha
[ptfret{1},~,counts{1},avgsig{1}] = portfolio_sort(hpr, signals_HF(:,:,1), opts);

% Skewness
[ptfret{2},~,counts{2},avgsig{2}] = portfolio_sort(hpr, signals_LF(:,:,2), opts);

% Skewness#2
[ptfret{3},~,counts{3},avgsig{3}] = portfolio_sort(hpr, signals_LF(:,:,3), opts);

% Bab
ptfret{4} = bab(hpr,signals_LF(:,:,4),rf);
ptfret{4} = bab(hpr,signals_HF(:,:,3),rf);

%% RV, RCOV, RBETA plots

% Alcoa
aa          = getTaqData('symbol','AA',[],[],[],'..\data\TAQ\sampled\5min\nobad_vw');
% issorted(aa.Datetime)
aa.Ret      = [NaN; diff(log(aa.Price))];
ion         = [true; diff(fix(aa.Datetime))~=0];
aa.Ret(ion) = NaN;

spy          = getSpy(5);
spy.Ret      = [NaN; diff(log(spy.Price))];
ion          = [true; diff(fix(spy.Datetime))~=0];
spy.Ret(ion) = NaN;

% Intersect dates
Date = intersect(fix(aa.Datetime), fix(spy.Datetime));
aa   = aa(ismember(fix(aa.Datetime), Date),:);
spy  = spy(ismember(fix(spy.Datetime), Date),:);

% NO OVERNIGHT
%%%%%%%%%%%%%%
rv   = @(subs) sqrt(accumarray(subs,spy.Ret.^2,[],@nansum));
rcov = @(subs,stock) accumarray(subs,spy.Ret.*stock.Ret,[],@nansum);

% Daily
[Date,~,iday] = unique(serial2yyyymmdd(spy.Datetime));
RVd           = rv(iday);
RCVd          = rcov(iday,aa);

% Monthly
[~,pos,imonth] = unique(serial2yyyymmdd(spy.Datetime)/100,'last');
RVm            = rv(imonth);
RCVm           = rcov(imonth,aa);

% Yearly, rolling every month
RVy      = NaN(size(RCVm));
RCVy     = NaN(size(RCVm));
rcov_run = @(idx,stock) nansum(spy.Ret(idx).*stock.Ret(idx));
for ii = 12:size(RCVm,1)
    idx      = ismember(imonth, ii-11:ii);
    RVy(ii)  = sqrt(nansum(spy.Ret(idx).^2));
    RCVy(ii) = rcov_run(idx,aa);
end

% WITH OVERNIGHT
%%%%%%%%%%%%%%%%
spy.Permno(:,1) = 84398;
spy             = addOvernightRet(spy);
aa              = addOvernightRet(aa);

rv   = @(subs) sqrt(accumarray(subs,spy.Ret.^2,[],@nansum));
rcov = @(subs,stock) accumarray(subs,spy.Ret.*stock.Ret,[],@nansum);

RVd_on  = rv(iday);
RCVd_on = rcov(iday,aa);

RVm_on  = rv(imonth);
RCVm_on = rcov(imonth,aa);

RVy_on   = NaN(size(RVm));
RCVy_on  = NaN(size(RVm));
rcov_run = @(idx,stock) nansum(spy.Ret(idx).*stock.Ret(idx));
for ii = 12:size(RVm,1)
    idx         = ismember(imonth, ii-11:ii);
    RVy_on(ii)  = sqrt(nansum(spy.Ret(idx).^2));
    RCVy_on(ii) = rcov_run(idx,aa);
end

dt_daily   = yyyymmdd2datetime(Date);
dt_monthly = serial2datetime(spy.Datetime(pos));

% Plot RV
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RVd,    dt_monthly, RVm,    dt_monthly, RVy,...
     dt_daily, RVd_on, dt_monthly, RVm_on, dt_monthly, RVy_on)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3,0.4],'YtickLabel',{'0','10\%','20\%','30\%','40\%'})
set(gca,'TickLabelInterpreter','latex')
print('SPYrv','-depsc','-r200')

% Plot RCV
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RCVd,    dt_monthly, RCVm,    dt_monthly, RCVy,...
     dt_daily, RCVd_on, dt_monthly, RCVm_on, dt_monthly, RCVy_on)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3],'YtickLabel',{'0','10\%','20\%','30\%'},'Ylim',[-0.01,0.3])
set(gca,'TickLabelInterpreter','latex')

% Plot RBETA
figure
set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.62],'PaperPositionMode','auto')
plot(dt_daily, RCVd./RVd.^2,    dt_monthly, RCVm./RVm.^2,    dt_monthly, RCVy./RVy.^2,...
     dt_daily, RCVd_on./RVd_on.^2, dt_monthly, RCVm_on./RVm_on.^2, dt_monthly, RCVy_on./RVy_on.^2)
legend('daily','monthly','yearly','d_on','m_on','y_on','Location','NorthWest')
legend boxoff
set(gca,'Ytick',[0,0.1,0.2,0.3],'YtickLabel',{'0','10\%','20\%','30\%'},'Ylim',[-0.01,0.3])
set(gca,'TickLabelInterpreter','latex')
