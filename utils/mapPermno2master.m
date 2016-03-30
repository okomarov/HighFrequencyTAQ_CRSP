function mst = mapPermno2master
% MAPPERMNO2MASTER Maps TAQ issues to CRSP's Permno

% Map CUSIP table with dates

%% 1. Import mst-cusip (TAQmaster)

mstcusip = loadresults('TAQmaster');
mstcusip = mstcusip(:,{'SYMBOL','CUSIP','FDATE'});

% 8-CUSIP
taqcusip = char(mstcusip.CUSIP);
taqcusip = taqcusip(:,1:8);

% Drop empty or null cusips
empty          = all(bsxfun(@eq, taqcusip, repmat(' ',1,8)),2);
null           = all(bsxfun(@eq, taqcusip, repmat('0',1,8)),2);
ikeep          = ~(empty | null);
mstcusip       = mstcusip(ikeep,:);
mstcusip.CUSIP = cellstr(taqcusip(ikeep,:));

% Unique records for mst-cusip
mstcusip = sortrows(mstcusip,{'CUSIP','FDATE'});
ikeep    = isfeatchange(mstcusip(:,{'CUSIP','SYMBOL','FDATE'}),3);
mstcusip = mstcusip(ikeep,:);

%% 2. Import CRSP data
msenames             = loadresults('msenames');
msenames             = msenames(msenames.Nameendt > 19923112,{'Permno','Namedt','Nameendt','Ncusip'});
idx                  = msenames.Namedt <= 19923112;
msenames.Namedt(idx) = 19930101;

% Drop empty
ikeep    = ~cellfun('isempty',msenames.Ncusip);
msenames = msenames(ikeep,:);

%% 3. Pre-filter CRSP and TAQ by cusip
idx      = ismember(msenames.Ncusip, mstcusip.CUSIP);
msenames = msenames(idx,:);
idx      = ismember(mstcusip.CUSIP, msenames.Ncusip);
mstcusip = mstcusip(idx,:);

%% 4. Import mst-symbol (mst)
master = load(fullfile('.\data\TAQ','master'),'-mat');
% Replace p with PR and remove . from TAQ symbols
ids    = master.ids;
ids    = regexprep(ids,'p','PR');
ids    = regexprep(ids,'\.','');
%% 5. Pre-filter mst-symbol and mst-cusip by symbol
idx       = ismember(ids, mstcusip.SYMBOL);
ids(~idx) = {''};
idx       = ismember(mstcusip.SYMBOL, ids);
mstcusip  = mstcusip(idx,:);

% Propagate symbol pre-filtering to CRSP
idx      = ismember(msenames.Ncusip, mstcusip.CUSIP);
msenames = msenames(idx,:);
%% 6. Map cusip to one set of numbers
msenames            = sortrows(msenames,{'Ncusip','Namedt'});
cusips              = intersect(msenames.Ncusip, mstcusip.CUSIP);
[~,msenames.Ncusip] = ismember(msenames.Ncusip, cusips);
[~,mstcusip.CUSIP]  = ismember(mstcusip.CUSIP, cusips);

%% 6. Cache (order matters)
% Cache mst-cusip in same order as ids
mstcusip               = sortrows(mstcusip,{'SYMBOL','FDATE'});
[unMstSymbol, ~, subs] = unique(mstcusip.SYMBOL);
nrows                  = accumarray(subs,1);
tmp                    = mat2cell(mstcusip, nrows);
[idx,pos]              = ismember(ids,unMstSymbol);
mstcusip               = cell(size(ids));
mstcusip(idx)          = tmp(pos(idx));

% Cache msenames by cusip
msenames = sortrows(msenames,{'Ncusip','Namedt'});
nrows    = accumarray(msenames.Ncusip,1);
msenames = mat2cell(msenames, nrows);

% Cache mst by Id
mst   = sortrows(master.mst,{'Id','Date'});
nrows = accumarray(mst.Id,1);
mst   = mat2cell(mst, nrows);

%% Map
cl     = class(msenames{1}.Permno);
Permno = arrayfun(@(x) zeros(x,1,cl),nrows,'un',0);
for ii = 1:numel(ids)
    disp(ii)
    symbol = ids{ii};
    if isempty(symbol)
        continue
    end
    
    % Fetch by cusip
    % NOTE: mst CUSIP is mapped into a progressive id and msenames is
    % already cached according to that number
    pos = unique(mstcusip{ii}.CUSIP,'stable');
    tmp = cat(1,msenames{pos});
    
    % Intersect CUSIP ranges 
    %
    % Example mstcusip{ii} and tmp:
    %
    %     SYMBOL    CUSIP     FDATE  
    %     ______    _____    ________
    %     'S'       19322    19930104
    %     'S'       20149    20050815
    %     'S'       20149    20051018
    %     
    %     Permno     Namedt     Nameendt    Ncusip
    %     ______    ________    ________    ______
    %     14322     19930101    20020101    19322
    %     14322     20020102    20040609    19322
    %     14322     20040610    20050324    19322
    %     39087     19930101    20010823    20149   out-of-range
    %     39087     20010824    20020101    20149   out-of-range
    %     39087     20020102    20040609    20149   out-of-range
    %     39087     20040610    20050814    20149   out-of-range
    %     39087     20050815    20051026    20149
    %     39087     20051027    20130710    20149
    edges            = [mstcusip{ii}.FDATE; 99999999];
    [~,~,bin]        = histcounts(tmp.Namedt,edges);
    izero            = bin == 0;
    [~,~,bin(izero)] = histcounts(tmp.Nameendt(izero),edges);
    tmp(bin == 0,:)  = [];
    bin(bin == 0,:)  = [];
    tmp              = tmp(mstcusip{ii}.CUSIP(bin) == tmp.Ncusip,:);
    
    % Match date bucket
    for r = 1:size(tmp,1);
        idx             = in(mst{ii}.Date, [tmp.Namedt(r), tmp.Nameendt(r)]);
        Permno{ii}(idx) = tmp.Permno(r);
    end
end
mst        = cat(1,mst{:});
mst.Permno = cat(1,Permno{:});

% Map back to order in master
[~,pos] = ismembIdDate(master.mst.Id,master.mst.Date,mst.Id,mst.Date);
mst     = mst(pos,{'Id','Permno','Date'});

%% Post-process duplicates 

% Duplicate permno - date
master.mst.Permno      = mst.Permno;
[un,~,subs]            = unique(master.mst(:,{'Permno','Date'}));
dup                    = un(accumarray(subs,1) > 1,:);
dup(dup.Permno == 0,:) = [];

% Duplicate records, add symbols
idx               = ismembIdDate(master.mst.Permno, master.mst.Date, dup.Permno, dup.Date);
master.mst        = master.mst(idx,:);
master.mst.Symbol = master.ids(master.mst.Id);
master            = sortrows(master.mst,{'Permno','Date','Symbol'});

% Mark id - symbol to drop (rule: symbol with suffix is dropped)
master.Idrop = false(size(master,1),1);
permnos      = unique(master.Permno);
% isolved      = false(size(permnos));

for ii = 1:numel(permnos)
    p           = permnos(ii);
    idx         = master.Permno == p;
    symbols     = unique(master.Symbol(idx));
    isSubString = all(strncmpi(symbols{1}, symbols(2:end), numel(symbols{1})));
    if isSubString
        idx          = ismember(master.Symbol, symbols(2:end));
        master.Idrop = master.Idrop | idx;
%         isolved(ii)  = true;
    else
        % Drop all if conflict not resolved
        % We could dwell into TAQ CUSIP alst 4 digits...not worth it only
        % for 130 permnos left
        idx          = ismember(master.Symbol, symbols);
        master.Idrop = master.Idrop | idx;
    end
end
% Drop permno mapping
idrop             = master.Idrop;
idrop             = ismembIdDate(mst.Id,mst.Date, master.Id(idrop), master.Date(idrop));
mst.Permno(idrop) = 0;

% % Duplicate permno - date
% master = load(fullfile('.\data\TAQ','master'),'-mat');
% master.mst.Permno      = mst.Permno;
% [un,~,subs]            = unique(master.mst(:,{'Permno','Date'}));
% dup                    = un(accumarray(subs,1) > 1,:);
% dup(dup.Permno == 0,:) = [];

save(fullfile('.\results\', sprintf('%s_%s.mat', datestr(now,'yyyymmdd_HHMM'),'masterPermno')), 'mst')