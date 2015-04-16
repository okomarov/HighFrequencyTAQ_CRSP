function scores = estimateCondAlpha(lookback, betas, rets, f, C)
% estimateCondAlpha(betas, rets)
% Sorts betas by Date and UnID

% Cond alphas
% From r_{i,t+1} = alpha_it + f_{t+1} * beta_it + e_{i,t+1}:
% 1) estimate the daily risk-adjusted return ^Z_{i,t+1} = r_{i,t+1} - f_{t+1} * ^beta_it
% 2) estimate the conditional alpha E(^Z_{i,t+1}|C_t)

% Unstack returns and betas
rets = rets(~isnan(rets.RetCC),{'Date','UnID','RetCC'});
rets = unstack(rets,'RetCC','UnID');
rets = sortrows(rets,'Date');

betas = betas(~isnan(betas.Beta),{'Date','UnID','Beta'});
betas = unstack(betas(:,{'Date','UnID','Beta'}),'Beta','UnID');
betas = sortrows(betas,'Date');

% Intersect returns and betas
[~,ia,ib] = intersect(getVariableNames(betas),getVariableNames(rets));
betas     = betas(:,ia);
rets      = rets(:,ib);

% Ensure f sorted by date
f = sortrows(f(:,{'Date','RetCC'}),'Date');

% Reference dates
refdates = intersect(rets.Date, f.Date);

% Time align
notrail = true;
% Sample at t-1
betas   = sampledates(betas,refdates-1,notrail);
C       = sampledates(C    ,refdates-1,notrail);
% Sample at t
rets    = sampledates(rets ,refdates  ,notrail);
f       = sampledates(f    ,refdates  ,notrail);

% Residual returns
nobs = size(rets,1);
Z    = NaN(nobs, width(rets));
for c = 2:width(rets)
    Z(:,c) = rets{:,c} - betas{:,c}.*f.RetCC;
end

% POINT 2)
scores = NaN(nobs, width(rets));
X      = C{:,2:end};
for r = lookback:nobs
    pos = r-lookback+1:r;
    x   = X(pos,:);
    for c = 2:width(rets)
        y = Z(pos,c);
        if nnz(~isnan(y)) > 10
            mdl = regstats(y,x,'linear','tstat');
            scores(r,c) = mdl.tstat.beta(1);
        end
    end
end

scores = array2table(scores,'VariableNames',getVariableNames(rets));
scores.Date = rets.Date;
end
