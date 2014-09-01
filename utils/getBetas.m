function [betas, unids] = getBetas(issp, iscs)
% GETBETAS Loads betas, applies sp500 and/or common shares filers, unstacks
%
%   NOTE: Does NOT fill in between NaNs

if nargin < 1 || isempty(issp), issp = false; end
if nargin < 2 || isempty(iscs), iscs = false; end

% Load all betas
loadresults('Betas','betas')
if isa(betas,'dataset'), betas = dataset2table(betas); end

% Filter for sp500 members
if issp
    idx   = issp500member(betas(:,{'Date','UnID'}));
    betas = betas(idx,:);
end

% Filter for common share (share type code 10 and 11)
if iscs
    idx   = iscommonshare(betas(:,{'UnID','Date'}));
    betas = betas(idx,:);
end

% Remove overlapping
[un,~,subs] = unique(betas(:,{'UnID','Date'}));
overlap     = un(accumarray(subs,1) > 1,:);
[~,idx]     = setdiff(betas(:,{'UnID','Date'}), overlap);
betas       = betas(idx,:);
% taq2crsp(ismember(taq2crsp.permno, taq2crsp.permno(taq2crsp.ID == overlap.ID(n))) | taq2crsp.ID == overlap.ID(n),:)

if nargout == 2
    unids = unique(betas.UnID);
end

% Unstack betas
betas = unstack(betas(:,{'Date','UnID','Beta'}), 'Beta','UnID');
betas = sortrows(betas,'Date');

% Convert to double
names = getVariableNames(betas);
betas = varfun(@double, betas);
betas = setVariableNames(betas, names);
end
