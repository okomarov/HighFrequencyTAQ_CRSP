function tb = getFFIndustryCodes(tb, numind, isPanel)
% GETFFINDUSTRYCODES Get Fama and French industry group codes for given Date - Permno pairs
%
%   GETFFINDUSTRYCODES(TB,NUMIND) TB is a table with Permno and Date in yyyymmdd format
%       and NUMIND specifies on how many industry codes to map.
%       Currently available classifications are 49 and 12 industries.

if nargin < 2
    error('getFFIndustryCodes:missingIndustry','Missing industry specification.')
end

if nargin < 3
    isPanel = false;
end

% Ensure it is sorted by id-date, i.e. some date changes correspond
% to same permno, we need to sort
idx = diff(tb.Date) == 0;
if ~all(tb.Permno(idx) ~= tb.Permno([false;idx]))
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

% Import mapping of SIC to the FF industry portfolios
switch numind
    case 49
        [sic2ffptf, desc] = getFF49Classification();
    case 12
        [sic2ffptf, desc] = getFF12Classification();
    otherwise
        error('getFFIndustryCodes:invalidIndustry','Invalid indusrty specification.')
end

% Map
for ii = 1:size(sic2ffptf,1)
    from            = sic2ffptf.Sic_from(ii);
    to              = sic2ffptf.Sic_to(ii);
    idx             = in(ccm.Sic, [from, to],'[]');
    ccm.FFid(idx,1) = uint8(sic2ffptf.FF_id(ii));
end

% Exclude unmapped
% 3990 - Referring to generic grouping?
% 6797 - Does not exist
% 9995 - NONCLASSIFIABLE ESTABLISHMENTS
% 9997 - Same as above (?)
idx        = ccm.FFid == 0;
ccm(idx,:) = [];

% Cache
N          = numel(unP_tb);
[~,pos]    = ismember(ccm.Lpermno, unP_tb);
ccm_cached = cache2cell(ccm{:,:}, pos,[],[N,1]);
[~,pos]    = ismember(tb.Permno, unP_tb);
tb_dates   = cache2cell(tb.Date, pos,[],[N,1]);

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