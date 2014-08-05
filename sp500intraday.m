loadresults('taq2crsp')
loadresults('spconst')
loadresults('dseshares')
path2data = '.\data\TAQ\sampled';
master = load(fullfile(path2data, 'master'), '-mat');

% Filter out non members
dseshares  = dseshares(ismember(dseshares.PERMNO, spconst.Id),:);
taq2crsp   = dataset2table(taq2crsp(ismember(taq2crsp.permno, spconst.Id),:));
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

% Permno to UnID
% NOTE: it's a 1 - many join 
dseshares = innerjoin(dseshares, taq2crsp, 'LeftKeys','Permno','RightKeys','permno',...
                                'RightVariables','ID');
% Shrink to matching in mst and take unique records
dseshares = unique(dseshares(ismember(dseshares.ID, master.mst.UnID),:));
                            
% Cache
dseshares = innerjoin(dseshares, master.mst, 'LeftKeys',{'ID','Date'},...
                                      'RightKeys',{'UnID','Date'},...
                                      'RightVariables','File');
                                  
dseshares = dseshares(:,{'Permno','Date','ID', 'Shrout','File'});
tickers = master.ids(unique(master.mst.Id));
data = getData(master, tickers, [],[],[],path2data); 

% clearvars -except dseshares path2data
% Analyze('spproxy',[],dseshares,fullfile(path2data,'S5m_*.mat'),1)

