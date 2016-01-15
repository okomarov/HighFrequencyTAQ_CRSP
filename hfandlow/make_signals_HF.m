function [signals,hpr,rf] = make_signals_HF(permno,date,factors)

reton = loadresults('return_intraday_overnight');

[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = size(data,2);

signals = NaN(nmonths, nseries,4);
hpr     = NaN(nmonths, nseries);
rf      = NaN(nmonths, 1);


for ii = 1:nmonths
    % Get prices
    from = unDt(ii)*100+1;
    to   = unDt(ii)*100+31;
    data = getTaqData('permno',permno,from,to,[],'..\data\TAQ\sampled\5min\nobad_vw');
    
    % Add returns
    data.Ret        = [NaN; diff(log(data.Price))];
    ion             = [true; diff(fix(data.Datetime)) ~= 0] |...
                      [true; diff(data.Permno) ~= 0];
    data.Ret(ion,1) = NaN;
    data            = addOvernightRet(data,reton);
    
    % Unstack
    data.Date   = serial2yyyymmdd(data.Datetime);
    data.Time   = serial2hhmmss(data.Datetime);
    data        = sortrows(unstack(data(:,{'Date','Time','Permno','Ret'}),'Ret','Permno'),{'Date','Time'});
    permnoFound = xstr2num(data.Properties.VariableNames(3:end));
    data        = data{:,3:end};

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

% Holding period return
for ii = 1:nmonths
    imonth = midx == ii;

    r                 = data(imonth,:);
    inan              = isnan(r);
    r(inan)           = 0;
    hpr(ii,:)         = prod(1+r)-1;
    hpr(ii,all(inan)) = NaN;

    rf(ii) = prod(1+factors.RF(imonth)/100)-1;
end
end