%% Options
OPT_LAG       = 1;
OPT_ALLVALUED = true;
%% Data 
dsf         = loadresults('dsf');
dsf.IsMicro = isMicrocap(dsf,'Prc',OPT_LAG);

ff = loadresults('F-F_Research_Data_5_Factors_2x3_daily_TXT');
ff = ff(ismember(ff.Date, unique(dsf.Date)),:);

ret    = sortrows(unstack(dsf(:,{'Permno','Date','Ret'}),'Ret','Permno'),'Date');
date   = ret.Date;
permno = ret.Properties.VariableNames(2:end);
ret    = ret{:,2:end};

signals = make_signals(ret,date,ff);