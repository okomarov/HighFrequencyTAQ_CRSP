%% Opts
OPT_RANGES = [930, 1000, 1030, 1100, 1130, 1200, 1230, 1300, 1330, 1400, 1430, 1500,1530]';
OPT_LAG    = 1;

OPT_STANDARDIZE = true;
OPT_VOL_AVG     = 'e';
OPT_LAG_VOL     = 40;
%% Data

% Returns
ret        = loadresults('halfHourRet930');
ret.T93000 = double(ret.T93000);
for ii = 2:numel(OPT_RANGES)
    tmp         = loadresults(sprintf('halfHourRet%d',OPT_RANGES(ii)));
    [~,pos]     = ismembIdDate(ret.Permno, ret.Date, tmp.Permno, tmp.Date);
    tname       = sprintf('T%d',OPT_RANGES(ii)*100);
    ret.(tname) = double(tmp.(tname)(pos));
end
ret = sortrows(ret,{'Permno','Date'});

% Variance
vol           = loadresults('halfHourVol930');
vol.RV5_93000 = double(vol.RV5_93000);
for ii = 2:numel(OPT_RANGES)
    tmp         = loadresults(sprintf('halfHourVol%d',OPT_RANGES(ii)));
    [~,pos]     = ismembIdDate(vol.Permno, vol.Date, tmp.Permno, tmp.Date);
    tname       = sprintf('RV5_%d',OPT_RANGES(ii)*100);
    vol.(tname) = double(tmp.(tname)(pos));
end
vol = sortrows(vol,{'Permno','Date'});

% isequal(vol.Permno, ret.Permno)
% isequal(vol.Date, ret.Date)

% Extract 
permnos = ret.Permno;
dates   = ret.Date;
ret     = ret{:,3:end};
vol     = vol{:,3:end};

% Group observations by dates, for clustered standard errors
[~,~,g] = unique(dates);
%% Sigma

% 40-day moving average lagged 1 standard deviation
sigma = sqrt(tsmovavg(vol,OPT_VOL_AVG, OPT_LAG_VOL,1));
sigma = [NaN(OPT_LAG,size(sigma,2)); sigma(1:end-OPT_LAG,:)]; % ex-ante

% Do not use lags from other permnos
SHIFT            = OPT_LAG_VOL-1+OPT_LAG;
igroup           = [false(SHIFT,1); permnos(1+SHIFT:end) == permnos(1:end-SHIFT)];
sigma(~igroup,:) = NaN;
%% Multivariate
% Regress
% h13 ~ constant + h11, h10 ... h1
% h12 ~ constant + h10, h19 ... h1
% ...

n              = numel(OPT_RANGES);
[tstat, coeff] = deal(NaN(n-1,n-2));
for ypos = n:-1:3
    % Standardizing returns by vol
    xpos = 1:ypos-2;
    if OPT_STANDARDIZE
        X = ret(:,xpos)./sigma(:,xpos);
        y = ret(:,ypos)./sigma(:,ypos);
    else
        X = ret(:,xpos);
        y = ret(:,ypos);
    end
    tb         = clusterreg(y,X,g,'linear');
    r          = 1:xpos(end)+1;
    c          = ypos-2;
    tstat(r,c) = tb.tStat;
    coeff(r,c) = tb.Estimate;
end
%% Univariate
% Regress
% h13 ~ constant + h11
% h13 ~ constant + h10
% ...
% h12 ~ constant + h10
% h12 ~ constant + h9
% ...

n              = numel(OPT_RANGES);
[tstat, coeff] = deal(NaN(n-1,n-2));
for ypos = n:-1:3
    for xpos = ypos-2:-1:1
        % Standardizing returns by vol
        if OPT_STANDARDIZE
            X = ret(:,xpos)./sigma(:,xpos);
            y = ret(:,ypos)./sigma(:,ypos);
        else
            X = ret(:,xpos);
            y = ret(:,ypos);
        end
        tb         = clusterreg(y,X,g,'linear');
        r          = xpos+1;
        c          = ypos-2;
        tstat(r,c) = tb.tStat(2);
        coeff(r,c) = tb.Estimate(2);
    end
end