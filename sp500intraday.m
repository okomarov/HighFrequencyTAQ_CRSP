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

% Data matrix of number of shares
sharesdata     = table2array(dseshares.Panel);

% Add permno to master
[~,pos]           = ismember(master.mst.UnID,taq2crsp.ID);
master.mst.Permno = taq2crsp.permno(pos);

% LOOP by date
[dates,~,subs] = unique(master.mst.Date);
ndates         = numel(dates);
res            = cell(ndates,1);
tic
pctRunOnAll warning off MATLAB:table:ModifiedVarnames
poolStartup(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
parfor ii = 1:ndates
    disp(ii)
    % Retrieve reference permnos for given date
    row    = ismembc2(dates(ii), spconst.Panel.Date);
    shares = sharesdata(row,2:end);
    
    % Get mst records
    idates = subs == ii;
    mst    = master.mst(idates,:);
    
    % Restrict to reference permnos (overlapping?)
    [~, imst, ishares] = intersect(mst.Permno, spconst.Id);
    mst = mst(imst,:);
    
    % Get data
    data = getTaqData(mst,[],[],[],[],path2data);
    
    % Use the backfill price strategy to build proxy (should update
    % previous day index with gradually changing prices, need overnight return)
    data.Price = flipud(nanfillts(data.Price(end:-1:1)));
    
    % Ensure all datetimes are to the minute
    dt = datevec(data.Datetime);
    data.Datetime = datenum(dt(:,1),dt(:,2),dt(:,3),dt(:,4),dt(:,5),0);
    
    % Pivot
    prices     = unstack(data(:,{'Permno','Datetime','Price'}),'Price','Permno');
    datetimes  = prices.Datetime;
    prices     = table2array(prices(:,2:end));
    
    % (should use pervious day close)
    capitaliz = double(shares(ishares)).* prices(1,:);
    index = sum(bsxfun(@times, prices, capitaliz./sum(capitaliz)),2);
    
    % Store results
    res{ii} = table(datetimes, index, 'VariableNames',{'Datetime','Price'});
end
pctRunOnAll warning on MATLAB:table:ModifiedVarnames
delete(gcp)
res = cat(1,res{:});

% Plot
loadresults('spysampled')

subplot(211)
plot(datetime(datevec(res.Datetime)), res.Price)
title 'sp500 proxy - rebalanced daily with open price'

subplot(212)
plot(datetime(datevec(spysampled.Datetime)), spysampled.Price)
title 'spyders'

saveas(gcf,'SP500proxy vs SPY.png')
toc