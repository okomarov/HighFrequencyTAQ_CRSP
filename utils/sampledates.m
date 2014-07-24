function tspanel = sampledates(tspanel, refdates)
% SAMPLEDATES Sample the time-series panel (unstacked table)
%
%   SAMPLEDATES(TSPANEL, REFDATES) 


% Union of dates
dates    = tspanel.Date;
alldates = union(dates, refdates);

% Map to union
[~,pos] = ismember(alldates,dates);

% Fill stretched periods with previous val
pos(pos == 0) = NaN;
pos           = nanfillts(pos);
tspanel       = tspanel(pos,:);
tspanel.Date  = alldates;

% Restrict to refdates
idx     = ismembc(alldates, refdates);
tspanel = tspanel(idx,:);
end