%% Options
OPT_VW = true;

OPT_LAG    = 1;
OPT_PTF_UN = 10;

%% Data 
dsf = loadresults('dsf');

myunstack = @(tb,vname) sortrows(unstack(tb(:,{'Permno','Date',vname}),vname,'Permno'),'Date');

% Returns
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

% Factors
ff = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
ff = ff(ismember(ff.Date, unique(dsf.Date)),:);
%% Signals
[signals, hpr, rf] = make_signals(ret,date,ff);
nsig           = size(signals,3);
snames         = {'alpha','skewh','skewr','betas'};
correlations   = corrxs(signals,snames);

%% Lag
% End-of-Month
signals = signals(1:end-OPT_LAG,:,:);
isMicro = isMicro(1:end-OPT_LAG,:);
cap     = cap(1:end-OPT_LAG,:);

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
[ptfret{1},~,counts{1},avgsig{1}] = portfolio_sort(hpr, signals(:,:,1), opts);

% Skewness
[ptfret{2},~,counts{2},avgsig{2}] = portfolio_sort(hpr, signals(:,:,2), opts);

% Skewness#2
[ptfret{3},~,counts{3},avgsig{3}] = portfolio_sort(hpr, signals(:,:,3), opts);

% Bab
ptfret{4} = bab(hpr,signals(:,:,4),rf);
