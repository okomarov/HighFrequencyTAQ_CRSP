function [betas, rets] = estimateCondAlpha(betas, rets)
% estimateCondAlpha(betas, rets)
% Sorts betas by Date and UnID

% Cond alphas
% From y_it = alpha_it + beta_it * f_t + e_it, estimate:
% 1) ^beta_it with time-series regression
% 2) ^f_t with cross-section regression
% 3) ^alpha_it = y_it - ^beta_it * ^f_t

% POINT 1)
% Get/estimate betas (input)
% Sort according to date
betas = sortrows(betas,{'Date','UnID'});

% POINT 2)
% Daily rets from return_overnight (input)
% Intersect sets, order guaranteed
[~,ia,ib] = intersect(rets(:,{'Date','UnID'}), betas(:,{'Date','UnID'}));
rets      = rets(ia,:);
betas     = betas(ib,:); 

% Subs by day 
[~,~,subs] = unique(rets.Date);

% Fhat
Fhat = CovOnVar(betas.Beta, rets.Totret, subs);

% POINT 3)
betas.Alpha = rets.Totret - Fhat(subs).*betas.Beta;
end

function beta = CovOnVar(x,y,subs)
% Low frequency
Exy  = accumarray(subs, x.*y, [],@nanmean);
Ex   = accumarray(subs,    x, [],@nanmean);
Ey   = accumarray(subs,    y, [],@nanmean);
Cov  = Exy - Ex.*Ey;
Var  = accumarray(subs,  x,[], @nanvar);
beta = Cov./Var;
end
