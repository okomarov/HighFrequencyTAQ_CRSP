function bool = isExchange(permnos, dates, code)

% ISEXCHANGE Checks from CRSP if and when stocks are traded on exchange(s)
% 
%   ISEXCHANGE(PERMNOS, DATES) PERMNOS can be 'xNumbers' or numeric and
%                              DATES should be in the yyyymmdd format.
%
%   ISEXCHANGE(..., CODE) CODE can be one or a combination from the list
%                         below. (default CODE = 1)
%
%   CRSP list of valid EXCHCD codes:
%
%         Code  Description
%          -2   Halted by Primary Listing Exchange
%          -1   Suspended by Primary Listing Exchange
%           0   Not Trading on Primary Listing Exchange
%           1   NYSE
%           2   NYSE MKT
%           3   NASDAQ
%           4   Arca
%           5   Mutual Funds (As Quoted By NASDAQ)
%          10   Boston Stock Exchange
%          13   Chicago Stock Exchange
%          16   Pacific Stock Exchange
%          17   Philadelphia Stock Exchange
%          19   Toronto Stock Exchange
%          20   Over-The-Counter (Non-NASDAQ Dealer Quotations)
%          31   When-Issued Trading on NYSE
%          32   When-Issued Trading on NYSE MKT
%          33   When-Issued Trading on NASDAQ
%          34   When-Issued Trading on Arca

%% Checks

if iscellstr(permnos)
    permnos = xstr2num(permnos);
end
if nargin < 3 || isempty(code)
    code = 1; 
else
    ivalid = ismember(code, cast([-2:5, 10,13,16,17,19,20,31:34],'like',code));
    if ~all(ivalid);
        error('isExchange:invalidCode','Invalid CODE. See help for a list of valid CRSP codes.')
    end
end

%% Engine
% Load CRSP info
mnames = loadresults('msenames');

% Filter out by exchange code
idx    = ismember(mnames.Exchcd, code);
mnames = mnames(idx, {'Permno','Namedt','Nameendt','Exchcd'});

% Filter out permnos
idx    = ismember(mnames.Permno, permnos);
mnames = mnames(idx,:);

% LOOP by batches of # permnos
N    = numel(permnos);
bool = false(numel(dates), N);
s    = 1000;
for ii = 1:s:N
    % Select relevant portion of permnos from msenames 
    range         = ii:min(ii+s-1,N);
    [idx,pos]     = ismember(mnames.Permno, permnos(range));
    range         = unique(range(pos(idx)));
    % Pivot from to date ranges of selected exchange codes and convert to logical
    exchcd        = pivotFromTo(mnames(idx,:), dates);
    bool(:,range) = logical(exchcd.Panel{:,2:end});
end
end