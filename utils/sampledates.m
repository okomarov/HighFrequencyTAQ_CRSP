function tspanel = sampledates(tspanel, refdates, notrail)
% SAMPLEDATES Sample the time-series panel (unstacked table)
%
%   SAMPLEDATES(TSPANEL, REFDATES, NOTRAIL) 

if nargin < 3, notrail = false; end
 
% if ~issorted(tspanel.Date)
%     error('sampledates:unsorted','Dates are not sorted in ascending order.')
% end

% Union of dates
dates    = tspanel.Date;
alldates = union(dates, refdates);

% Map to union
[~,pos] = ismember(alldates,dates);

% Fill stretched periods with previous val
nullpos      = pos == 0;
pos(nullpos) = NaN;
pos          = nanfillts(pos);
tspanel      = tspanel(pos,:);
tspanel.Date = alldates;
if notrail
    from = find(~nullpos,1,'last')+1;
    tspanel(from:end,2:end) = {NaN};
end

% Restrict to refdates
idx     = ismember(alldates, refdates);
tspanel = tspanel(idx,:);
end