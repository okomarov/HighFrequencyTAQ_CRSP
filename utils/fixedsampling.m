function [newprices, newdates, actdates] = fixedsampling(dates,prices, grid)

% FIXEDSAMPLING Sample time series at specific points in time
%
%   ... = FIXEDSAMPLING(DATES,PRICES,GRID)
%
%       DATES   double vector of dates:
%
%       PRICES  numeric vector of prices. Must have the same number of
%               elements as DATES.
%
%       GRID    a vector of intraday times, i.e. values in [0 1], sorted in
%               ascending order. The grid will be replicated for all days in
%               DATES.
%               Note: 1 second = 1/86400. 18:39 = (3600* 18+39 *60)/86400.
%
%   [NEWPRICES, NEWDATES, ACTDATES] = ...
%
%       NEWPRICES   sampled prices at gridpoints with last-price interpolation
%       NEWDATES    gridpoint dates
%       ACTDATES    actual dates picked for each gridpoint
%
% Examples:
%
%  See also:

% Author: Oleg Komarov (oleg.komarov@hotmail.it)
% Tested on R2013a.
% 13 Dec 2013 - Created
%  4 Jul 2014 - Changed in/out arguments, simplified, cleaned up, added warnings


% Ninput
narginchk(3,3)

% Dates
ndates = numel(dates);
if ndates ~= numel(prices) || ~isvector(dates) || ~isvector(prices)
    error('fixedsampling:data','DATES and PRICES should be vectors with equal number of elements.')
end
if isa(dates,'single')
    warning('fixedsampling:singleDates','Single DATES have at most 1.5 hour precision, use double.')
end
if isempty(prices)
    error('fixedsampling:emptyPrices','PRICES are empty.')
end
if ~issorted(dates)
    error('fixedsampling:sortedDates','DATES should be sorted in ascending order.')
end

% Interval
if ~isnumeric(grid) && ~issorted(grid) && any(grid < 0 | grid > 1)
    error('fixedsampling:interval','GRID should be a vector of values in [0 1].')
end
if isrow(grid)
    grid = grid';
end

% -------------------------------------------------------------------------
% ENGINE: create fixed grid
% -------------------------------------------------------------------------

% Slightly less than millisecond tolerance
tol = 1/(60*60*24*1000)-eps;
% Position of last observations for each day
last    = [0; find(diff(fix(dates))); ndates];
nseries = numel(last)-1;

% Add beginning of day to grid to avoid sample overflow to next day
% NOTE: should work even if 0 is already there
grid = [zeros(1,class(grid)); grid];

% Vector of intraday times [0 1]: replicate grid for every day
ngrid = numel(grid);
grid  = bsxfun(@plus,grid,fix(dates(last(2:end),1)).');

% Subtract tolerance to keep beginning of day in next bucket
grid(1,:) = grid(1,:)-tol;
grid      = grid(:);

% Bin with lb < x <= ub
n = histc(dates, grid + tol);

% -------------------------------------------------------------------------
% ENGINE: select observations from the grid with last-price interpolation
% -------------------------------------------------------------------------

% Preallocate
newprices      = NaN(ngrid,nseries);
% Which bins gridpoints have new data, i.e. after previous gridpoint
idx            = logical(n);
% Position of new data, skipping initial 0s
pos            = cumsum(n(idx));
newprices(idx) = prices(pos);
% Carry over data from last available gridpoint
newprices      = nanfillts(newprices);
% Skip very firs gridpoint (artificially placed there) and exclude very bin which includes all point >= ub(end)
newprices      = reshape(newprices(1:end-1),ngrid*nseries-1,1);
% Exclude beginning-of-day to first gridpoint artificial overflow bucket
newprices(ngrid:ngrid:end,:) = [];

% New dates
if nargout >= 2
    grid(1:ngrid:end,:) = [];
    newdates = grid;
end

% Actual dates
if nargout == 3
    actdates      = NaN(ngrid,nseries);
    actdates(idx) = dates(pos);
    actdates      = reshape(nanfillts(actdates(1:end-1)),ngrid*nseries-1,1);
    actdates(ngrid:ngrid:end,:) = [];
end
end