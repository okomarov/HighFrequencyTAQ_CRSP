function tb = getFF49IndustryCodes(tb, isPanel)
% GETFF49INDUSTRYCODES Get Fama and French 49 industry group codes for given Date - Permno pairs
%  
%   GETFF49INDUSTRYCODES(TB) TB is a table with Permno and Date in yyyymmdd format

if nargin < 2, isPanel = false; end

% Ensure it is sorted by id-date, i.e. some date changes correspond
% to same permno, we need to sort
idx = diff(tb.Date) == 0;
if ~all(tb.Permno(idx) ~= tb.Permno([false;idx]));
    [tb,isort] = sortrows(tb,{'Permno','Date'});
    SORT_BACK  = true;
else
    SORT_BACK = false;
end

% Compustat sic codes
try
    ccm = loadresults('ccmsic');
catch
    ccm = loadresults('ccmsic','..\results');
end
ccm = sortrows(ccm(:,{'Lpermno','Linkdt','Linkenddt','Sic'}),1:2);
ccm = consolidateFromTo(ccm);

% Intersect permnos
unP_tb = unique(tb.Permno);
idx    = ismember(ccm.Lpermno, unP_tb);
ccm    = ccm(idx,:);

% Import mapping of SIC to the 49 FF industry portfolios
[sic2ffptf, desc] = getClassificationFF49();

% Map
for ii = 1:size(sic2ffptf,1)
    from            = sic2ffptf.Sic_from(ii);
    to              = sic2ffptf.Sic_to(ii);
    idx             = in(ccm.Sic, [from, to],'[]');
    ccm.FFid(idx,1) = uint8(sic2ffptf.Id_FF49(ii));
end

% Exclude unmapped
% 3990 - Referring to generic grouping?
% 6797 - Does not exist
% 9995 - NONCLASSIFIABLE ESTABLISHMENTS 
% 9997 - Same as above (?)
idx        = ccm.FFid == 0;
ccm(idx,:) = [];

% Cache
N               = numel(unP_tb);
ccm_cached      = cell(1,numel(unP_tb));
ccm             = ccm{:,:};
[~,pos]         = ismember(unique(ccm(:,1)), unP_tb);
ccm_cached(pos) = cache2cell(ccm, ccm(:,1));
tb_dates        = cache2cell(tb.Date,tb.Permno);

FFid = cell(N,1);
for ii = 1:N
    ccm_ith      = ccm_cached{ii};
    tb_dates_ith = tb_dates{ii};
    FFid{ii}     = zeros(size(tb_dates_ith,1),1,'uint8');
    
    if isempty(ccm_ith)
        continue
    end
    
    % Assign ff classification
    for r = 1:size(ccm_ith,1)
        idx           = in(tb_dates_ith, ccm_ith(r,2:3));
        FFid{ii}(idx) = ccm_ith(r,end);
    end
    
    % Backfill as in Goyal 2014 - No Size Anomalies in U.S. Bank Stock Returns
    izero = FFid{ii} == 0;
    if any(izero) && ~all(izero)
        tmp        = double(FFid{ii});
        tmp(izero) = NaN;
        FFid{ii}   = int8(flipud(nanfillts(tmp(end:-1:1))));
    end
end
tb.FFid = cat(1,FFid{:});

% Unstack
if isPanel
    tb = sortrows(unstack(tb(:,{'FFid','Permno','Date'}),'FFid','Permno'),'Date');
elseif SORT_BACK
    tb(isort,:) = tb;
end
end