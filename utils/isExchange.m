function bool = isExchange(data, code)

% ISEXCHANGE Checks from CRSP if and when stocks are traded on exchange(s)
% 
%   ISEXCHANGE(DATA) DATA should be structure with the exact fields:
%                       .Panel  - with numeric data and NaNs
%                       .Permno - a cell array of 'xNumbers' (strings)
%                       .Date   - in yyyymmdd numeric format
%
%                    DATA can also be a fints or a table panel.
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

if nargin < 2 || isempty(code)
    code = 1; 
else
    ivalid = ismember(code, cast([-2:5, 10,13,16,17,19,20,31:34],'like',code));
    if ~all(ivalid);
        error('isExchange:invalidCode','Invalid CODE. See help for a list of valid CRSP codes.')
    end
end

% DATA handling
if isa(data, 'fints')
    tmp  = fts2mat(data,1);
    data = struct('Panel' , tmp(:,2:end),...
                  'Permno', {fieldnames(data,1)},...
                  'Date',   serial2yyyymmdd(tmp(:,1)));
elseif isa(data,'table')
    vnames = getVariableNames(data);
    data   = struct('Panel',data{:,2:end},'Permno',{vnames(2:end)},'Date',data.Date);
elseif isfield(data,'Data')
    data.Panel = data.Data;
    data       = rmfield(data,'Data');
end
% Convert xPermnos
data.Permno = xstr2num(data.Permno);

% Load CRSP info
mnames = loadresults('msenames');

% Filter out by exchange code
idx    = ismember(mnames.EXCHCD, code);
mnames = mnames(idx, {'PERMNO','NAMEDT','NAMEENDT','EXCHCD'});

% LOOP by batches of 5000 permnos
bool = false(size(data.Panel));
s    = 1000;
N    = numel(data.Permno);
for ii = 1:s:N
    % Select relevant portion of permnos from msenames 
    range       = ii:min(ii+s-1,N);
    [idx,pos]   = ismember(data.Permno(range), mnames.PERMNO);
    pos         = pos(idx);
    tmpnames    = mnames(pos,:);
    % Pivot from to date ranges of selected exchange codes and convert to logical 
    exchcd      = pivotFromTo(tmpnames, data.Date);
    bool(:,idx) = logical(exchcd.Panel{:,2:end});
end
end