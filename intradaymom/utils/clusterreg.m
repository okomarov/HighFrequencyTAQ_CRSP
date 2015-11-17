function [res, varBhat] = clusterreg(y, X, g, model)
% CLUSTERREG OLS regression with one-way or two-way clustered standard errors
%
%   CLUSTERREG(y, X, g)
%       where y is the dependent variable, X is the matrix of regressors,
%       and g is a vector of group indices. All inputs must have same
%       number of observations (rows) and g can have two columns with a
%       second group of clustering indices
%
%   RES = clusterreg(...)
%       RES is a table with:
%           .Betas
%           .Se         clustered standard errors 
%           .Tstats     t-statistics
%
%   References: 
%   [1] Cameron, A. C., J. B. Gelbach, and D. L. Miller. "Robust inference 
%       with multiway clustering." Journal of Business & Economic 
%       Statistics vol.29, no. 2 (2011)

if nargin < 4, model = 'linear'; end
X = x2fx(X,model);

nonnan = ~any(isnan(X),2);
if any(~nonnan)
    X = X(nonnan,:);
    y = y(nonnan,:);
    g = g(nonnan,:);
end

Betas = X\y;
e     = y - X*Betas;

% Cluster robust variance on first group of indices
varBhat = clusteredVar(X, e, g(:,1));

if size(g,2) == 2
    % With two clustering indices, the final variance is equal to the sum
    % of the variances clustered by each group, minus the variance from the
    % intersected clusters, i.e. B_hat = VCV_g1 + VCV_g2 - VCV_g12.
    varBhat = varBhat + clusteredVar(X, e, g(:,2)) - clusteredVar(X, e, g);
end

% Calculate standard errors and t-stats
Se     = sqrt(diag(varBhat));
Tstats = Betas ./ Se;

% Return the calculated values
res = table(Betas, Se, Tstats);
end

function varBhat = clusteredVar(X, e, g)
[N, k] = size(X);

if size(g,2) == 2
    [G,~,glabel] = unique(g, 'rows');
else
    [G,~,glabel] = unique(g);
end
M = numel(G);

X_g = cache2cell(X,glabel);
e_g = cache2cell(e,glabel);
B   = 0;
for ii = 1:M
    disp(ii)
    B = B + (X_g{ii}'*e_g{ii})*(e_g{ii}'*X_g{ii});
end

% Calculate cluster-robust variance matrix estimate
q_c     = (N-1)/(N-k)*M/(M-1);
XX      = X'*X;
varBhat = q_c * (XX\B)/XX; % inv(X'*X)*B*inv(X'*X)
end