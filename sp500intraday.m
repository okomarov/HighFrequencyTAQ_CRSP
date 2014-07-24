loadresults('taq2crsp')
loadresults('spconst')
loadresults('dseshares')
path2data = '.\data\TAQ\sampled';
load(fullfile(path2data, 'master'), '-mat');

% Filter out non members
dseshares = dseshares(ismember(dseshares.PERMNO, spconst.Id),:);
taq2crsp  = taq2crsp(ismember(taq2crsp.permno, spconst.Id),:);

% Time consolidation
idx       = isfeatchange(dseshares(:,{'PERMNO', 'SHROUT','SHRSDT'}));
from      = find(idx);
to        = [from(2:end)-1; numel(idx)];
dseshares = [dseshares(from,{'PERMNO','SHROUT','SHRSDT'}), dseshares(to,'SHRENDDT')];
dseshares = pivotFromTo(dseshares(:,{'PERMNO','SHRSDT','SHRENDDT','SHROUT'}));

% Sample at reference dates
refdates        = unique(mst.Date);
dseshares.Panel = sampledates(dseshares.Panel,refdates);
spconst.Panel   = sampledates(spconst.Panel,refdates);

% Filter out by membership
vnames = getVariableNames(spconst.Panel);
for c = 2:size(spconst.Panel,2)
    field = vnames{c};
    dseshares.Panel.(field)(~spconst.Panel.(field)) = 0;
end

% Stack dseshares back to tall dataset
permnos   = repmat(dseshares.Id, size(dseshares.Panel,1),1);
dseshares = stack(dseshares.Panel,vnames(2:end),...
                  'NewDataVariableName','Shrout');
dseshares.Permno = permnos;
dseshares = dseshares(dseshares.Shrout ~= 0, {'Date','Permno','Shrout'});


% Just start building index by looping each day and exclude based on the
% ticker

% Substitute permno in dseshares with unId, then cache it with file number
% and write routine in analyze

[idx, pos] = ismember(mst.UnID,taq2crsp.ID);
Pivot([taq2crsp.permno(pos(idx)), mst.Date(idx), ones(nnz(idx),1)],[],1,0);
