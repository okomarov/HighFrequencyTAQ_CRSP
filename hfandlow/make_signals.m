function [signals,hpr] = make_signals(ret,date,factors)
[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = size(ret,2);

signals = NaN(nmonths, nseries,4);
hpr     = NaN(nmonths, nseries);
for ii = 1:nmonths
    imonth = midx == ii;
    nobs   = nnz(imonth);

    % Month of all non-nan returns
    r      = ret(imonth,:);
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

for ii = 12:nmonths
    iyear = ismember(midx, ii-12+1:ii);
    nobs  = nnz(iyear);

    % Month of all non-nan returns
    r      = ret(iyear,:);
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
    imonth            = midx == ii;
    r                 = ret(imonth,:);
    inan              = isnan(r);
    r(inan)           = 0;
    hpr(ii,:)         = prod(1+r)-1;
    hpr(ii,all(inan)) = NaN;
end
end