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
