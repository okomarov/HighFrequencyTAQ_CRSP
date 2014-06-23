function [out,actdates] = fixedsampling(data,scheme,grid)

% FIXEDSAMPLING Sample time series at specific points in time
%
%   ... = FIXEDSAMPLING(DATA,SCHEME,GRID)
%
%       DATA    is an m by 2 double matrix:
%                   - column 1      increasing serial dates
%                   - column 2      values
%
%
%       SCHEME  is a string indicating how the datapoints corresponding to
%               the gridpoints are selected (we are always in Calendar Time).
%
%               Only available option for now is:
%                   - 'Previous'    last available observation up to the
%                                 actual gridpoint (aka last price
%                                 interpolation).
%
%                Eventually available:
%                   - 'First'       first observation after the previous
%                                 gridpoint.
%                   - 'Last'        last observation after the previous
%                                 gridpoint.
%                   - 'Min'         observation with minimum price in the
%                                 interval (previous actual] gridpoint.
%                   - 'Max'         observation with maximum price in the
%                                 interval (previous actual] gridpoint.
%                   - 'Next'        first available observation after the
%                                 actual gridpoint (aka first price
%                                 interpolation).
%                   - 'Linear'      linear interpolation between datapoints
%                                 selected with 'Previous' and 'Next'.
%                   - 'Nearest'     closest datapoint to the gridpoint.
%                   - 'Uniform'     daily first and last observations are
%                                 included by default and the remaining
%                                 INTERVAL-2 points are selected as
%                                 'Nearest' from a uniformly spread grid.
%
%
%       GRID    a vector of intraday times, i.e. values in [0 1], sorted in
%               ascending order. The grid will be replicated for all days in
%               DATA.
%               Note: 1 second = 1/86400. 18:39 = (3600* 18+39 *60)/86400.
%
%   [OUT,ACTDATES] = ...
%
%       OUT         filtered DATA with fixed grid times.
%       ACTDATES    actual dates picked for each gridpoint
%
% Examples:
%
%  See also:

% Author: Oleg Komarov (oleg.komarov@hotmail.it)
% Tested on R2013a.
% 13 Dec 2013 - Created

% Ninput
narginchk(3,3)

% Data
szData = size(data);
if isempty(data) || ~isfloat(data) || szData(2) ~= 2
    error('fltprice:data','DATA should be a single/double m by 2 matrix.')
end
if ~issorted(data(:,1))
    error('fltprice:data1stColumn','DATA''s 1st column (serial dates) should be sorted in ascending order.')
end

% Scheme
if ~ischar(scheme) || ~isvector(scheme)
    error('fltprice:scheme','SCHEME should be a string.')
else
    % Try to match
    [err,scheme] = getScheme(scheme);
    if ~isempty(err)
        error('fltprice:scheme',err)
    end
end

% Interval
if ~isnumeric(grid) && ~issorted(grid) && any(grid < 0 | grid > 1)
    error('fltprice:interval','GRID should be a vector of values in [0 1].')
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
last    = [0; find(diff(fix(data(:,1)))); szData(1)];
nseries = numel(last)-1;

% Add beginning of day to grid to avoid sample overflow to next day
% NOTE: should work even if 0 is already there
grid = [zeros(1,class(grid)); grid];

% Vector of intraday times [0 1]: replicate grid for every day
ngrid = numel(grid);
grid  = bsxfun(@plus,grid,fix(data(last(2:end),1)).');

% Subtract tolerance to keep beginning of day in next bucket
grid(1,:) = grid(1,:)-tol;
grid      = grid(:);

% Bin with lb < x <= ub
n = histc(data(:,1), grid + tol);

% -------------------------------------------------------------------------
% ENGINE: select observations from the grid according to the SCHEME
% -------------------------------------------------------------------------

switch scheme
    %     case 'First'
    %         % If no values are found within a grid interval we don't want to
    %         % interpolate the first value. This is accomplished with consequent
    %         % indexing of the same position instead of using position themselves.
    %         idx = false(szData(1),1);
    %         n   = cumsum([2; n(1:end-1)]);
    %         % Take care values that fall to the next day
    %         n(ismember(n,last+1) | n > szData(1)) = [];
    %         idx(n) = true;
    %
    %     case 'Last'
    %         % Same concept for missing values as in First
    %         idx = false(szData(1),1);
    %         idx(cumsum([1; n(1:end-1)])) = true;
    %
    %     case 'Linear'
    %         % Calculate Previous and Next and interpolate time and price
    %         prev = cumsum([1; n(1:end-1)]);
    %         next = prev + 1;
    %         % Remove first and last values
    %         prev(ismember(prev,last+1)) = [];
    %         [idx,loc] = ismember(prev,last);
    %         prev = prev(~idx);
    %         next(ismember(next,[last+1;last+2]) | next > szData(1)) = [];
    %         % Output
    %         out = (data(prev,:) + data(next,:))/2;
    %         % Add first/last back and sort back (not efficient but will do for
    %         % now)
    %         out = sort([out;
    %                     data(last(1:end-1)+1,:);
    %                     data(last(setdiff(loc,0)),:)]);
    %
    %     case 'Max'
    %         if ~strcmp(type,'FixedTime')
    %             [n,bin] = histc(data(:,1), grid + tol);
    %         end
    %         % Only prices because there can be multiple max prices per grid
    %         % interval
    %         out = accumarray([1;bin(2:end-1)],data(1:sum(n)+1,2),[],@max);
    %
    %     case 'Min'
    %         if ~strcmp(type,'FixedTime')
    %             [n,bin] = histc(data(:,1), grid + tol);
    %         end
    %         % Only prices because there can be multiple min prices per grid
    %         % interval
    %         out = accumarray([1;bin(2:end-1)],data(1:sum(n)+1,2),[],@min);
    %     case 'Nearest'
    %         % Previous and next
    %         n = cumsum([1; n(1:end-1)]);
    %         n = [n, n+1];
    %         % If next overshots the day set to last observation
    %         if n(end,2) > szData(1)
    %            n(end,2) = szData(1);
    %         end
    %         % Find minimum distance from gridpoint (nearest)
    %         [~,pos] = min(diff([data(n(:,1),1) grid(:)...
    %                             data(n(:,2),1)],[],2),[],2);
    %         numn = size(n,1);
    %         % Select previous or next whichever closer
    %         idx  = n((1:numn).' + (pos-1)*numn);
    %
    %     case 'Next'
    %         % First price interpolation
    %         % Use positions directly to obtain carry-on interpolation
    %         idx = cumsum([2; n(1:end-1)]);
    %         idx(ismembc(idx,last+1) | idx > szData(1)) = [];
    %
    %     case 'Standard'
    %         idx = mcolon(last(1:end-1)+1, last(2:end), grid);
    
    % Last price interpolation
    case 'Previous'
        % Preallocate
        [tmp{1:max(nargout,1)}] = deal(NaN(ngrid,nseries));
        % Which bins gridpoints have new data, i.e. after previous gridpoint
        idx         = logical(n);
        % Position of new data, skipping initial 0s
        pos         = cumsum(n(idx));
        tmp{1}(idx) = data(pos,2);
        % Carry over data from last available gridpoint
        tmp{1}      = nanfillts(tmp{1});
        % Skip very firs gridpoint (artificially placed there) and exclude very bin which includes all point >= ub(end)
        out         = [grid(2:end), reshape(tmp{1}(1:end-1),ngrid*nseries-1,1)];
        if nargout == 2
            tmp{2}(idx)           = data(pos,1);
            actdates              = reshape(nanfillts(tmp{2}),ngrid*nseries,1);
            actdates(ngrid:ngrid:end) = [];
        end
        % Exclude beginning-of-day to first gridpoint artificial overflow bucket
        out(ngrid:ngrid:end,:) = [];
end
end

% getScheme ---------------------------------------------------------------
function [err,scheme] = getScheme(scheme)
% Initialize error
err = '';
% Available schemes
whichScheme = {'First','Last','Linear','Max','Min','Nearest','Next','Previous','Standard'};
% Try to match
idx  = strncmpi(scheme,whichScheme,numel(scheme));
% # of matches
nidx = nnz(idx);

% Ambiguous scheme
if  nidx == 2
    err = sprintf('SCHEME ''%s'' is ambiguous. Did you mean ''%s'' or ''%s''?',scheme,whichScheme{idx});
    % No match
elseif nidx == 0
    err = sprintf('SCHEME ''%s'' unrecognized.',scheme);
    % Regular match
else
    scheme = whichScheme{idx};
end
end