function [signals,hpr,rf] = make_signals_HF(permno,date,factors)

reton = loadresults('return_intraday_overnight');

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = numel(permno);

signals = NaN(nmonths, nseries,4);
hpr     = NaN(nmonths, nseries);
rf      = NaN(nmonths, 1);

% Get market return
mkt             = getSpy(5, min(date),max(date));
mkt.Permno(:,1) = 84398;
mkt             = price2ret(mkt);
mkt             = addOvernightRet(mkt,reton);

% Get rf


for ii = 1:nmonths
    [ret, permnoFound] = getHighFreqRet(unDt(ii),permno,reton);

    % Month of all non-nan returns
    r      = data(imonth,:);
    nonans = all(~isnan(r));
    ngood  = nnz(nonans);
    r      = r(:,nonans);

    % Alpha
    Y     = bsxfun(@minus,r, factors.RF(imonth)/100);
    X     = [ones(nobs,1), factors.MktMinusRF(imonth)/100];
    coeff = NaN(2,ngood);

    for jj = 1:ngood
        coeff(:,jj) = X\Y(:,jj);
    end
    signals(ii,nonans,1) = coeff(1,:);

    % Skewness
    signals(ii,nonans,2) = skewness(r);
    signals(ii,nonans,3) = sqrt(nobs) * sum(r.^3) ./ sum(r.*r).^1.5;
end

% Betas
for ii = 12:nmonths
    iyear = ismember(midx, ii-12+1:ii);
    nobs  = nnz(iyear);

    % Month of all non-nan returns
    r      = data(iyear,:);
    inan   = isnan(r);
    nonans = sum(~inan) >= 200;
    ngood  = nnz(nonans);
    inan   = inan(:,nonans);
    r      = r(:,nonans);

    % Betas
    Y     = bsxfun(@minus,r, factors.RF(iyear)/100);
    X     = [ones(nobs,1), factors.MktMinusRF(iyear)/100];
    coeff = NaN(2,ngood);

    for jj = 1:ngood
        idx         = ~inan(:,jj);
        coeff(:,jj) = X(idx,:)\Y(idx,jj);
    end
    signals(ii,nonans,4) = coeff(2,:);
end
end

function [data, permnoFound] = getHighFreqRet(month, permno, reton)
    % Get prices
    from = month*100+1;
    to   = month*100+31;
    data = getTaqData('permno',permno,from,to,[],'..\data\TAQ\sampled\5min\nobad_vw');

    % Add returns
    data = price2ret(data);
    data = addOvernightRet(data,reton);

    % Unstack
    data.Date   = serial2yyyymmdd(data.Datetime);
    data.Time   = serial2hhmmss(data.Datetime);
    data        = sortrows(unstack(data(:,{'Date','Time','Permno','Ret'}),'Ret','Permno'),{'Date','Time'});
    permnoFound = xstr2num(data.Properties.VariableNames(3:end));
    data        = data{:,3:end};
end

function data = price2ret(data)
    data.Ret        = [NaN; diff(log(data.Price))];
    ion             = [true; diff(fix(data.Datetime)) ~= 0] |...
                      [true; diff(data.Permno) ~= 0];
    data.Ret(ion,1) = NaN;
end