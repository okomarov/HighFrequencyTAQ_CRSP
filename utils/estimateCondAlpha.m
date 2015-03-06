function alpha = estimateCondAlpha(lookback, freq, useovern,  useproxy, sp500only, commononly)
% estimateCondAlpha(lookback, freq, useovern,  useproxy, sp500only, commononly)

% Cond alphas
% From y_it = alpha_it + beta_it * f_t + e_it, estimate:
% 1) ^beta_it with time-series regression
% 2) ^f_t with cross-section regression
% 3) ^alpha_it = y_it - ^beta_it*^f_t

% POINT 1)
% Get/estimate betas
beta = getBetas(lookback, freq, useovern, ~useproxy, sp500only, commononly, true);
% Sort according to date
beta = sortrows(beta,{'Date','UnID'});

% POINT 2)
% Daily rets [TODO use the dsfquery from crsp]
rets = loadresults('dailyret');

% Intersect sets, order guaranteed
[~,ia,ib] = intersect(rets(:,{'Date','UnID'}), beta(:,{'Date','UnID'}));
rets      = rets(ia,:);
beta      = beta(ib,:); 

% Subs by day
[~,~,subs] = unique(rets.Date);

% Numerator
Exy2 = accumarray(subs, (beta.Beta .* rets.Dret).^2,[],@nansum);
Ex   = accumarray(subs, beta.Beta ,[],@nanmean);
Ey   = accumarray(subs, rets.Dret    ,[],@nanmean);
Cov  = Exy2 - Ex.*Ey;
% Denominator
Ex2  = accumarray(subs, beta.Beta.*2 ,[],@nansum);
% Fhat 
Fhat = Cov./Ex2;

% POINT 3)
beta.Alpha = rets.Dret - Fhat(subs).*beta.Beta;
alpha = beta;
end