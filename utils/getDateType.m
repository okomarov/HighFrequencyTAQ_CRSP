function dttype = getDateType(date)
% GETDATETYPE Parse date input and return its type
%
%   DATE can be:
%       date               - 'scalar' numeric date in yyyymmdd format
%       inf (or [inf,inf]) - 'all' dates
%       [inf date]         - 'le', i.e. takes all dates up to a threshold (<= date)
%       [date inf]         - 'ge', takes all dates after a threshold i.e. >= date
%       [date  date]       - 'fromto', i.e. takes all dates in between
%       [date; date; ...]  - 'set', i.e. only the specified dates
%
% See also: isYYYYMMDD

sz   = size(date);
iinf = isinf(date);

if ~(isYYYYMMDD(date) | iinf)
    error('getDateType:invalidDate','DATE should be in YYYYMMDD format.')
end

% [inf, inf] case
if isrow(date) && sz(2) == 2 && all(iinf)
    date = inf;
end

if isscalar(date)

    % Inf case
    if iinf
        dttype = 'all';
        % Single date
    else
        dttype = 'scalar';
    end

% [from, to]
elseif sz(2) == 2 && all(~iinf)
    dttype = 'fromto';

% [inf, to] or [from, inf] cases
elseif sz(2) == 2 && any(iinf)
    if iinf(1)
        dttype = 'le';
    elseif iinf(2)
        dttype = 'ge';
    end

% Column of dates
elseif sz(2) == 1
    if any(iinf)
        error('A set of dates cannot have INFs.')
    end
    dttype = 'set';
else
    error('DATE cannot have more than 2 columns.')
end
end