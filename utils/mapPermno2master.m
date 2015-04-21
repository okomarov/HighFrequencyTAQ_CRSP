function mst = mapPermno2master
% Import CRSP data
msenames             = loadresults('msenames');
msenames             = msenames(msenames.NAMEENDT > 19923112,{'PERMNO','NAMEDT','NAMEENDT','NCUSIP','TICKER','TSYMBOL'});
idx                  = msenames.NAMEDT <= 19923112;
msenames.NAMEDT(idx) = 19930101;
% Import TAQ data
taqmaster            = loadresults('TAQmaster');
taqmaster            = taqmaster(:,{'SYMBOL','CUSIP','FDATE'});

% Msenames cusip
msecusip = msenames.NCUSIP;

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
% [unNcusip,~,msenames.NCUSIP] = unique(msenames.NCUSIP);
% pivotFromTo(msenames(:,{'PERMNO','NAMEDT','NAMEENDT','NCUSIP'}));

% % Consolidate TAQ (from dates)
% taqmaster = sortrows(taqmaster(:,{'CUSIP','SYMBOL','FDATE'}),{'CUSIP','FDATE'});
% idx       = isfeatchange(taqmaster,3);
% taqmaster = taqmaster(idx,:);

% Consolidate MSE (from dates)
msenames          = sortrows(msenames,{'PERMNO','NAMEDT'});
pstart            = find(isfeatchange(msenames(:,{'PERMNO','TSYMBOL','NCUSIP','NAMEDT'}),[1,4]));
pend              = [pstart(2:end)-1; size(msenames,1)];
To                = msenames.NAMEENDT(pend);
msenames          = msenames(pstart,:);
msenames.NAMEENDT = To;

% Fill in empty symbols with ticker
idx                   = strcmpi(msenames.TSYMBOL,'');
msenames.TSYMBOL(idx) = msenames.TICKER(idx);

% Forward fill intermediate and end of period empty symbols
msenames     = sortrows(msenames,{'NCUSIP','NAMEDT'});
iempty       = strcmpi(msenames.TSYMBOL,'');
pos          = NaN(size(iempty));
pos(~iempty) = find(~iempty);
pos          = nanfillts(pos);

% Backward fill beginning of period empty symbols
[~,pstart]       = unique(msenames.PERMNO);
pstart           = intersect(pstart,find(iempty));
pos(pstart)      = pstart + 1;
msenames.TSYMBOL = msenames.TSYMBOL(pos);

% Load master
master = load(fullfile('.\data\TAQ','master'),'-mat');

% Replace p with PR and remove . from TAQ symbols
ids = master.ids;
ids = regexprep(ids,'p','PR');
ids = regexprep(ids,'\.','');

% Intersect symbols
symbols          = intersect(ids, msenames.TSYMBOL);
idx              = ismember(ids,symbols);
master.ids(~idx) = {[]};
idx              = ismember(msenames.TSYMBOL, symbols);
msenames         = sortrows(msenames(idx,:),{'TSYMBOL','NAMEDT'});

% Cache mst
mst   = sortrows(master.mst,{'Id','Date'});
nrows = accumarray(mst.Id,1);
mst   = mat2cell(mst, nrows);

% Preallocate 
Permno = cell(numel(nrows),1);

for ii = 1:numel(master.ids)
    symbol     = master.ids{ii};
    Permno{ii} = zeros(nrows(ii),1,'like',msenames.PERMNO);
    if isempty(symbol)
        continue
    end
    
    % CRSP symbol table 
    isymbol = strcmpi(symbol, msenames.TSYMBOL);
    tmp     = msenames(isymbol, :);
        
    % Match date bucket
    [~,ptmp] = histc(mst{ii}.Date, [tmp.NAMEDT; 99999999]);
    imatched = ptmp ~= 0;
    ptmp     = ptmp(imatched);
    
    % Assign permno
    Permno{ii}(imatched) = tmp.PERMNO(ptmp);
end
mst        = cat(1,mst{:});
mst.Permno = cat(1,Permno{:});

% Map back to order in master
[~,pos] = ismembIdDate(mst.Id,mst.Date,master.mst.Id,master.mst.Date);
mst     = mst(pos,{'Id','Permno','Date'});

save(fullfile('.\results\', sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),'masterPermno')), 'mst')