function [betas, rets] = estimateCondAlpha(betas, rets, f)
% estimateCondAlpha(betas, rets)
% Sorts betas by Date and UnID

% Cond alphas
% From r_{i,t+1} = alpha_it + f_{t+1} * beta_it + e_{i,t+1}:
% 1) estimate the daily risk-adjusted return ^Z_{i,t+1} = r_{i,t+1} - f_{t+1} * ^beta_it
% 2) estimate the conditional alpha E(^Z_{i,t+1}|C_t)

% Unstack returns and betas
rets = unstack(rets(:,{'Date','UnID','RetCC'}),'RetCC','UnID');
rets = sortrows(rets,'Date');

betas = unstack(betas(:,{'Date','UnID','Beta'}),'Beta','UnID');
betas = sortrows(betas,'Date');

% Intersect returns and betas
[~,ia,ib] = intersect(getVariableNames(betas),getVariableNames(rets));
betas = betas(:,ia);
rets = rets(:,ib);

% Ensure f sorted by date
f     = sortrows(f    ,'Date');

% Point betas to next day returns
[idx,pos] = ismember(betas.Date, f.Date(1:end-1));
from      = find(idx, 1,'first')-1;
betas     = betas(from:end,:);
pos       = pos(from:end);
f         = f.RetCC(pos+1);

[~  ,pos] = ismember(betas.Date, rets.Date(1:end-1,:));
rets      = rets(pos+1,:);

% Residual returns
Z = rets;
for ii = 2:width(rets)
    Z{:,ii} = rets{:,ii} - betas{:,ii}.*f;
end


% POINT 2)

betas.Alpha = rets.Totret - Fhat(subs).*betas.Beta;
end