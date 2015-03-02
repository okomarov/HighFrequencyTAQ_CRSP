function [betas, unids] = getBetas(lookback, freq, useon, useproxy, issp, iscs, keeplong)
% GETBETAS Loads betas, applies sp500 and/or common shares filers, unstacks
%
%   getBetas(lookback, freq, useon, useproxy, issp, iscs, keeplong)
%
%   NOTE: Does NOT fill in between NaNs

if nargin < 1 || isempty(lookback), lookback = 1;       end
if nargin < 2 || isempty(freq),     freq     = 5;       end
if nargin < 3 || isempty(useon),    useon    = true;    end
if nargin < 4 || isempty(useproxy), useproxy = false;   end
if nargin < 5 || isempty(issp),     issp     = false;   end
if nargin < 6 || isempty(iscs),     iscs     = false;   end
if nargin < 7 || isempty(keeplong), keeplong = false;   end

% Load all betas
betas = estimateBetas(lookback, freq, useon, useproxy);

% Filter for sp500 members
if issp
    fprintf('%s: filtering for sp500 members.\n', mfilename)
    idx   = issp500member(betas(:,{'Date','UnID'}));
    betas = betas(idx,:);
end

% Filter for common share (share type code 10 and 11)
if iscs
    fprintf('%s: filtering for common shares only.\n', mfilename)
    idx   = iscommonshare(betas(:,{'UnID','Date'}));
    betas = betas(idx,:);
end

% Remove overlapping
[un,~,subs] = unique(betas(:,{'UnID','Date'}));
overlap     = un(accumarray(subs,1) > 1,:);
if ~isempty(overlap)
    [~,idx] = setdiff(betas(:,{'UnID','Date'}), overlap);
    betas   = betas(idx,:);
end
% taq2crsp(ismember(taq2crsp.permno, taq2crsp.permno(taq2crsp.ID == overlap.ID(n))) | taq2crsp.ID == overlap.ID(n),:)

if nargout == 2
    unids = unique(betas.UnID);
end

% Unstack betas
if ~keeplong
    betas = unstack(betas(:,{'Date','UnID','Beta'}), 'Beta','UnID');
    betas = sortrows(betas,'Date');
end
    
% Convert to double
names = getVariableNames(betas);
betas = varfun(@double, betas);
betas = setVariableNames(betas, names);
end