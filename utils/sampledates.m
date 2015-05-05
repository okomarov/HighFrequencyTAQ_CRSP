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

% Expand
nullpos      = pos == 0;
pos(nullpos) = NaN;
pos          = nanfillts(pos);
from         = find(~isnan(pos),1,'first');
tspanel      = tspanel(pos(from:end),:);
if ~isempty(from) && from > 1
    filler  = NaN(from-1, width(tspanel)); 
    filler  = array2table(filler,'VariableNames', getVariableNames(tspanel)); 
    tspanel = [filler; tspanel];
end
tspanel.Date = alldates;

% Fill values
tspanel(:,2:end) = tbextend.varfun(@(x) nanfillts(x,notrail), tspanel(:,2:end),...
                   'RenameVariables',false);

% Eventually cut off the notrail part
if notrail
    from = find(~nullpos,1,'last')+1;
    if ~isempty(from) && from <= size(tspanel,1)
        tspanel(from:end,2:end) = {NaN};
    end
end

% Restrict to refdates
idx     = ismember(alldates, refdates);
tspanel = tspanel(idx,:);
end