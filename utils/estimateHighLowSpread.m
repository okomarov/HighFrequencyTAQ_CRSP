function hl = estimateHighLowSpread(LAG)
% ESTIMATEHIGHLOWSPREAD Estimates intraday spread with daily high and low prices
%
%   ESTIMATEHIGHLOWSPREAD(LAG) Minimum number of past days for the
%          time-series expectation
%
%   The formula for the spread is the following (see reference):
%     S     = 2(e^alpha - 1)/(1 + e^alpha)
%     alpha = [sqrt(2*beta) - sqrt(beta)]/[3 - 2*sqrt(2)] - sqrt(gamma/[3 - 2*sqrt(2)])
%     beta  = E[sum_j(ln(high_j/low_j)^2 for j = 0:1)] and expectation as average over past LAG days
%     gamma = [ln(high(t,t+1)/low(t,t+1))]^2 where high and low are rolling over two days
%
% Reference: 2012 Corwin, Schults - A Simple Way to Estimate Bid-Ask Spreads from Daily High and Low Prices - JF

% Load high low prices and cache
try
    hl = loadresults('estimateHighLowPrice');
catch
    hl = loadresults('estimateHighLowPrice','..\results');
end
hl      = hl(hl.Permno ~= 0,:);
invalid = isnan(hl.High) | isnan(hl.Low) | hl.High == 0 | hl.Low == 0;
hl      = hl(~invalid,:);

% Adjust by overnight
try
    reton = loadresults('return_intraday_overnight','hfandlow\results');
catch
    reton = loadresults('return_intraday_overnight','..\hfandlow\results');
end
[~,pos]  = ismembIdDate(hl.Permno, hl.Date, reton.Permno, reton.Date);
hl.HighA = double(hl.High) .* (1+reton.RetCO(pos));
hl.LowA  = double(hl.Low) .* (1+reton.RetCO(pos));
clear reton

gvar = findgroups(hl.Permno);
rows = accumarray(gvar,1);
hl   = cache2cell(hl, findgroups(hl.Permno));

for ii = 1:numel(hl)
    if rows(ii) < LAG
        [hl{ii}.Beta, hl{ii}.Gamma, hl{ii}.Alpha, hl{ii}.Spread] = deal(NaN(rows(ii),1));
    else
        hl{ii}.Beta        = cumsum(log(hl{ii}.HighA ./ hl{ii}.LowA).^2,1); % Moving sum
        hl{ii}.Beta(3:end) = hl{ii}.Beta(3:end) - hl{ii}.Beta(1:end-2);
        hl{ii}.Beta(1)     = NaN;
        hl{ii}.Beta(2:end) = movmean(hl{ii}.Beta(2:end),[LAG-1,0]);         % Increasing expectations
        hl{ii}.Beta(1:LAG) = NaN;

        twoday_high  = max(hl{ii}.HighA(1:end-1), hl{ii}.HighA(2:end));
        twoday_low   = min(hl{ii}.LowA(1:end-1), hl{ii}.LowA(2:end));
        hl{ii}.Gamma = [NaN; log(twoday_high./twoday_low).^2];

        den          = 3 - 2*sqrt(2);
        hl{ii}.Alpha = (sqrt(2*hl{ii}.Beta) - sqrt(hl{ii}.Beta))./den - sqrt(hl{ii}.Gamma./den);

        hl{ii}.Spread = 2*(exp(hl{ii}.Alpha)-1)./(1+exp(hl{ii}.Alpha));
    end
end
hl = cat(1, hl{:});
hl = hl(:,{'Permno','Date','Spread'});
end
