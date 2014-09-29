function mktcap = estimateMktCap(UnIDs,refdates)
% Load data
taq2crsp = loadresults('taq2crsp');
dseshares = loadresults('dseshares');
dseshares = unique(dseshares(:,{'PERMNO', 'SHRSDT','SHRENDDT','SHROUT'}));
path2data = '.\data\TAQ\sampled';
master    = load(fullfile(path2data, 'master'), '-mat');

% No permno, no shares
taq2crsp  = taq2crsp(~isnan(taq2crsp.permno),:);
imember   = ismember(taq2crsp.ID, UnIDs);
taq2crsp  = taq2crsp(imember,:);
UnIDs     = intersect(taq2crsp.ID, UnIDs);
permnos   = uint32(unique(taq2crsp.permno));
dseshares = dseshares(ismember(dseshares.PERMNO, permnos),:);

% Filter out non members from mst
master.mst = master.mst(ismember(master.mst.UnID, UnIDs),:);

% Add permno to mst (need it for getTaqData)
[~,pos] = ismember(master.mst.UnID, taq2crsp.ID);
master.mst.Permno = taq2crsp.permno(pos);

% Time consolidation
idx       = isfeatchange(dseshares(:,{'PERMNO', 'SHROUT','SHRSDT'}));
from      = find(idx);
to        = [from(2:end)-1; numel(idx)];
dseshares = [dseshares(from,{'PERMNO','SHROUT','SHRSDT'}), dseshares(to,'SHRENDDT')];

% Pivoting
dseshares = pivotFromTo(dseshares(:,{'PERMNO','SHRSDT','SHRENDDT','SHROUT'}));

% Time sampling
dseshares.Panel = sampledates(dseshares.Panel,refdates,true);

% Cache
master = accumarray(master.mst.File,(1:size(master.mst))',[],@(x) {master.mst(x,:)});
% Retrieve end of day prices
poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
N = numel(master);
price = cell(N,1);
parfor f = 1:N
    disp(f)
    price{f} = getTaqData(master{f},[],[],[],'Price', path2data);
    price{f} = price{f}(diff(price{f}.Datetime) < 0,:);
end
delete(gcp)
price = cat(1,price{:});
price.Date = uint32(serial2yyyymmdd(price.Datetime));
price.Datetime = [];

% Unstack price
unid2permno = unique(price(:,{'UnID','Permno'}));
unid2permno.Properties.RowNames = matlab.lang.makeValidName(cellstr(num2str(unid2permno.UnID)));
unid2permno.Permno = matlab.lang.makeValidName(cellstr(num2str(unid2permno.Permno)));
price = unstack(price(:,{'UnID','Date','Price'}),'Price','UnID');
price = sortrows(price, 'Date');

save debugstate

% Time sampling
price = sampledates(price,refdates,true);

% LOOP by UnID to create weights
for ii = 2:size(price,2)
    unid = price.Properties.VariableNames{ii};
    permno = unid2permno(unid,'Permno').Permno;
    shares = double(dseshares.Panel.(permno{:}));
    shares(shares == 0) = NaN;
    price.(unid) = shares .* double(price.(unid)([1 1:end-1]));
end
mktcap = price;
end