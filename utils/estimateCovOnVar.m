function beta = estimateCovOnVar(y, x, subs, isLowFreq)
% In a linear regression y = a + b*x + e
% beta = cov(x,y)/var(x)

% Implement lookback so that with high freq is simpoly a running sum of num
% and denom while for low frequency is lisghtly more complex

if nargin < 3, subs      = ones(size(x)); end
if nargin < 4, isLowFreq = false;         end

if isLowFreq
    % LF
    Exy2 = accumarray(subs, (x.*y).^2, [],@nansum);
    Ex   = accumarray(subs,         x, [],@nanmean);
    Ey   = accumarray(subs,         y, [],@nanmean);
    Cov  = Exy2 - Ex.*Ey;
    Var = accumarray(subs,  x,[], @nanvar);
else
    % HF
    Cov = accumarray(subs, x.*y, [],@nansum);
    Var = accumarray(subs, x.*2, [],@nansum);
end

beta = Cov./Var;
end
