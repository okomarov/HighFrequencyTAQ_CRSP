function mst = mapPermno2master
% Import CRSP data
msenames             = loadresults('msenames');
msenames             = msenames(msenames.Nameendt > 19923112,{'Permno','Namedt','Nameendt','Ncusip','Ticker','Tsymbol'});
idx                  = msenames.Namedt <= 19923112;
msenames.Namedt(idx) = 19930101;

% Import TAQ data
taqmaster = loadresults('TAQmaster');
taqmaster = taqmaster(:,{'SYMBOL','CUSIP','FDATE'});

% Msenames cusip
msecusip = msenames.Ncusip;

% Extract TAQ 8-cusip
taqcusip = char(taqmaster.CUSIP);
taqcusip = taqcusip(:,1:8);

% Find null and empty cusips
emptytaq = find(all(bsxfun(@eq, taqcusip, repmat(' ',1,8)),2));
emptymse = find(cellfun('isempty',msecusip));
null     = find(all(bsxfun(@eq, taqcusip, repmat('0',1,8)),2));
taqcusip = cellstr(taqcusip);

% Match cusips
[~,ia,ib] = intersect(msecusip, taqcusip);

% Exclude empty and null
ia = setdiff(ia, emptymse);
ib = setdiff(ib, union(emptytaq,null));

% Filter out
msenames        = msenames(ismember(msecusip, msecusip(ia)),:);
[idx,pos]       = ismember(taqcusip, taqcusip(ib));
taqmaster       = taqmaster(idx,:);
taqcusip        = taqcusip(ib);
taqmaster.CUSIP = taqcusip(pos(idx));

% % Check multiple ncusips for a permno (OK - no repeated dates)
% [unNcusip,~,msenames.Ncusip] = unique(msenames.Ncusip);
% pivotFromTo(msenames(:,{'Permno','Namedt','Nameendt','Ncusip'}));

% % Consolidate TAQ (from dates)
% taqmaster = sortrows(taqmaster(:,{'Cusip','Symbol','Fdate'}),{'Cusip','Fdate'});
% idx       = isfeatchange(taqmaster,3);
% taqmaster = taqmaster(idx,:);

% Consolidate MSE (from dates)
msenames          = sortrows(msenames,{'Permno','Namedt'});
pstart            = find(isfeatchange(msenames(:,{'Permno','Tsymbol','Ncusip','Namedt'}),[1,4]));
pend              = [pstart(2:end)-1; size(msenames,1)];
To                = msenames.Nameendt(pend);
msenames          = msenames(pstart,:);
msenames.Nameendt = To;

% Fill in empty symbols with ticker
idx                   = strcmpi(msenames.Tsymbol,'');
msenames.Tsymbol(idx) = msenames.Ticker(idx);

% Forward fill intermediate and end of period empty symbols
msenames     = sortrows(msenames,{'Ncusip','Namedt'});
iempty       = strcmpi(msenames.Tsymbol,'');
pos          = NaN(size(iempty));
pos(~iempty) = find(~iempty);
pos          = nanfillts(pos);

% Backward fill beginning of period empty symbols
[~,pstart]       = unique(msenames.Permno);
pstart           = intersect(pstart,find(iempty));
pos(pstart)      = pstart + 1;
msenames.Tsymbol = msenames.Tsymbol(pos);

% Load master
master = load(fullfile('.\data\TAQ','master'),'-mat');

% Replace p with PR and remove . from TAQ symbols
ids = master.ids;
ids = regexprep(ids,'p','PR');
ids = regexprep(ids,'\.','');

% Intersect symbols with ismember (beware of duplicates after previous step)
idx       = ismember(ids, msenames.Tsymbol);
ids(~idx) = {''};
idx       = ismember(msenames.Tsymbol, ids);
msenames  = sortrows(msenames(idx,:),{'Tsymbol','Namedt'});

% Cache mst by Id
mst   = sortrows(master.mst,{'Id','Date'});
nrows = accumarray(mst.Id,1);
mst   = mat2cell(mst, nrows);

% Preallocate 
Permno = cell(numel(nrows),1);

for ii = 1:numel(master.ids)
    symbol     = ids{ii};
    Permno{ii} = zeros(nrows(ii),1,'like',msenames.Permno);
    if isempty(symbol)
        continue
    end
    
    % CRSP symbol table 
    isymbol = strcmpi(symbol, msenames.Tsymbol) | strcmpi(symbol, msenames.Ticker);
    tmp     = msenames(isymbol, :);
        
    % Match date bucket
    [~,ptmp] = histc(mst{ii}.Date, [tmp.Namedt; 99999999]);
    imatched = ptmp ~= 0;
    ptmp     = ptmp(imatched);
    
    % Assign permno
    Permno{ii}(imatched) = tmp.Permno(ptmp);
end
mst        = cat(1,mst{:});
mst.Permno = cat(1,Permno{:});

% Map back to order in master
[~,pos] = ismembIdDate(mst.Id,mst.Date,master.mst.Id,master.mst.Date);
mst     = mst(pos,{'Id','Permno','Date'});

save(fullfile('.\results\', sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),'masterPermno')), 'mst')