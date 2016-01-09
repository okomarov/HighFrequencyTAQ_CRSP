function signals = make_signals(ret,date,factors)
[unDt,~,midx] = unique(date/100);
nmonths       = numel(unDt);
nseries       = size(ret,2);

signals = NaN(nmonths, nseries,3);
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
end