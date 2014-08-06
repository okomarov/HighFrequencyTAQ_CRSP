loadresults('taq2crsp')
loadresults('spconst')
loadresults('dseshares')
path2data = '.\data\TAQ\sampled';
master = load(fullfile(path2data, 'master'), '-mat');

% Filter out non members 
dseshares  = dseshares(ismember(dseshares.PERMNO, spconst.Id),:);
taq2crsp   = dataset2table(taq2crsp(ismember(taq2crsp.permno, spconst.Id),:));

% Drop tickers with suffixes
[~,~,subs] = unique(taq2crsp.permno);
N          = numel(subs);
idrop      = false(N,1);

for ii = 1:N
    % Index permno and corresponding tickers
    isubs   = subs == ii;
    symbols = taq2crsp.symbol(isubs);
    symlen  = taq2crsp.symlen(isubs);
    nsymb   = nnz(isubs);
        
    if nsymb > 1
        idx = false(nsymb, 1);
        for jj = 1:nsymb
            if all(idx), break,end
            % Mark to drop if substring
            idx = idx | ~cellfun('isempty', regexp(symbols, sprintf('^%s\\w+',symbols{jj}),'once'));
        end
        if any(idx), idrop(isubs) = idx; end
    end
end
taq2crsp(idrop,:) = [];
master.mst = master.mst(ismember(master.mst.UnID, taq2crsp.ID),:);

% Time consolidation
idx       = isfeatchange(dseshares(:,{'PERMNO', 'SHROUT','SHRSDT'}));
from      = find(idx);
to        = [from(2:end)-1; numel(idx)];
dseshares = [dseshares(from,{'PERMNO','SHROUT','SHRSDT'}), dseshares(to,'SHRENDDT')];
dseshares = pivotFromTo(dseshares(:,{'PERMNO','SHRSDT','SHRENDDT','SHROUT'}));

% Sample at reference dates
refdates        = unique(master.mst.Date);
dseshares.Panel = sampledates(dseshares.Panel,refdates);
spconst.Panel   = sampledates(spconst.Panel,refdates);

% Filter out by membership
vnames = getVariableNames(spconst.Panel);
for c = 2:size(spconst.Panel,2)
    field = vnames{c};
    dseshares.Panel.(field)(~spconst.Panel.(field)) = 0;
end

% Stack dseshares back to tall dataset (use nshares as implicit membership)
permnos   = repmat(dseshares.Id, size(dseshares.Panel,1),1);
dseshares = stack(dseshares.Panel,vnames(2:end),...
                  'NewDataVariableName','Shrout');
dseshares.Permno = permnos;
dseshares = dseshares(dseshares.Shrout ~= 0, {'Date','Permno','Shrout'});

% Get price data
tickers = master.ids(unique(master.mst.Id));
data    = getTaqData(master, tickers, [],[],[],path2data); 

% Map back data to permno
data.Id

% Check overlapping
[un,~,subs] = unique([data.Id fix(data.Datetime)]);
overlap     = un(accumarray(subs,1) > 1,:);
[~,idx]     = setdiff(Betasd(:,1:2), overlap);

% Unstack
data    = unstack(data, 'Price','Id');
