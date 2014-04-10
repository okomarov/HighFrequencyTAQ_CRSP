function tspanel = nanfillts(tspanel)

% NANFILLTS Fill consecutive NaNs in the time series with last available observation
%
%   NANFILLTS(TSPANEL) TSPANEL should be a numeric matrix (or a vector) 
%                      where each series is a column with the oldest 
%                      observations in the first row (element) and the 
%                      most recent ones in the last row (element).
%
%
%   TSPANEL = NANFILLTS(...) Returns the same matrix with filled NaNs.
%
%   NOTE: leading NaNs are left untouched.
%
% See also

% Author: Oleg Komarov (o.komarov11@bcsprime.com)
% Tested on R2013a Win7 64bit
% 13 Dec 2013 - Created

%% Checks
narginchk(1,1)

% Terminate if not applicable
if isempty(tspanel) || (isscalar(tspanel) && isnumeric(tspanel))
    return
end

% TSPANEL
sz = size(tspanel);
if ~isnumeric(tspanel) || numel(sz) > 2
    error('prevprice:numTspanel','TSPANEL should be a numeric vector or matrix.')
end
isrowvector = sz(1) == 1;
if isrowvector 
    tspanel = tspanel';
    sz      = sz([2,1]); 
end

%% Engine

% Pad with a row of NaNs at the beginning and end
bread    = NaN(1,sz(2));
sandwich = [bread; tspanel; bread];

% Find starting ending positions of consecutive NaNs in price. The double padding ensures that consecutive
% NaNs and the beginning and at the end also have a start/end pair
bread = false(1,sz(2));
inan  = [bread; isnan(tspanel); bread];
sten  = diff(inan);

% Start end row/col subs. I need row/col subs granularity instead of plain positions because I need to take
% into account the previous padding
[rst,~  ] = find(sten ==  1);
[ren,cen] = find(sten == -1);

% Preallocate row subs mapping to sandwich
row         = double(~inan);
% New row length with padding
lrow        = sz(1)+2;
% Convert the ends of each NaN sequence to positions, taking into account the padding and diff effect
endpos      = (cen-1)*lrow + ren+1;
% Place the lengths of NaN sequences + increment to next row
row(endpos) = ren-rst+1;
% Get row map to sandwich
row         = cumsum(row);

% Convert row map to position in sandwich
pos      = bsxfun(@plus,(0:sz(2)-1)*lrow,row);
sandwich = sandwich(pos);

% Extract only the forward-filled-NaN price matrix (i.e. exclude padding)
tspanel = sandwich(2:end-1,:);

% Eventually convert back to row vector
if isrowvector 
    tspanel = tspanel'; 
end
end