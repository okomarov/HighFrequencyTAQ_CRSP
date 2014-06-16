function tspanel = nanfillts(tspanel, notrail)

% NANFILLTS Fill NaNs in the time series with the previous valid observation
%
%   NANFILLTS(TSPANEL) TSPANEL should be a numeric matrix (or a vector) 
%                      where each series is a column with the oldest 
%                      observations in the first row (element) and the 
%                      most recent ones in the last row (element).
%
%   NANFILLTS(..., NOTRAIL) To avoid extrapolation of the last value of a 
%                           time series till the end of the panel, set 
%                           NOTRAIL to true (1). By default is false (0).
%
%
%   TSPANEL = NANFILLTS(...) Returns the same matrix with filled NaNs.
%
%   NOTE: leading NaNs are left untouched.
%
%
% Example:
%
%   % Input data       % Simple call       % No trailing
%                      nanfillts(ts)       nanfillts(ts,1)
%    ts =              ans =               ans = 
%        10  NaN            10  NaN             10  NaN
%       NaN  NaN            10  NaN             10  NaN
%       NaN    1            10    1             10    1
%         9    2             9    2              9    2
%       NaN  NaN             9    2            NaN  NaN
%
% Additional features:
% - <a href="matlab: web('https://github.com/okomarov/nanfillts','-browser')">GITHUB/nanfillts page</a>
%
% See also: INTERP1Q

% Author: Oleg Komarov (o.komarov11@imperial.ac.uk)
% Tested on R2013a Win7 64bit
% 13 Dec 2013 - Created
% 16 Jun 2014 - Improved memory footprint; added NOTRAIL; expanded help

%% Checks
narginchk(1,2)

% Terminate if not applicable
if isempty(tspanel) || (isscalar(tspanel) && isnumeric(tspanel))
    return
end

% TSPANEL
sz = size(tspanel);
if ~isnumeric(tspanel) || numel(sz) > 2
    error('nanfillts:numTspanel','TSPANEL should be a numeric vector or matrix.')
end
isrowvector = sz(1) == 1;
if isrowvector 
    tspanel = tspanel';
    sz      = sz([2,1]); 
end

if nargin < 2 || isempty(notrail)
    notrail = false;
elseif ~isscalar(notrail) || ~(notrail == 1 || notrail == 0 || islogical(notrail))
    error('nanfillts:invalidNotrail','NOTRAIL should be true (1) or false (0).')
end

%% Engine

% Find starting ending positions of consecutive NaNs in price. The double padding ensures that consecutive
% NaNs and the beginning and at the end also have a start/end pair
bread = false(1,sz(2));
inan  = [bread; isnan(tspanel); bread];
sten  = diff(int8(inan));

% Start end row/col subs. I need row/col subs granularity instead of plain positions because I need to take
% into account the previous padding
[rst,~  ] = find(sten ==  1);
[ren,cen] = find(sten == -1);

% Preallocate row subs mapping to sandwich and cast to smaller type
[rst,ren,cen,row,sz] = cast2smaller(rst,ren,cen,~inan,sz);

% Eventually don't trail last available value if time series ends earlier
if notrail
    idx      = ren > sz(1);
    ren(idx) = rst(idx);
end

% New row length with padding
lrow          = sz(1)+2;
% Convert the ends of each NaN sequence to positions, taking into account the padding and diff effect
endpos        = (cen-1)*lrow + ren+1;
% Place the lengths of NaN sequences + increment to next row
row(endpos)   = ren-rst+1;
% Get row map to sandwich and reduce indexing by 1
row           = cumsum(row)-1;
% Restore rows pointing to 0 (previous slice of bread)
row(row == 0) = 1;
row           = row(2:end-1,:);

% Convert row map to position in tspanel
pos     = bsxfun(@plus,(0:sz(2)-1)*sz(1),row);
% Extract only the forward-filled-NaN price matrix (i.e. exclude padding)
tspanel = tspanel(pos);

% Eventually convert back to row vector
if isrowvector 
    tspanel = tspanel'; 
end
end

function [rst,ren,cen,nonnan,sz] = cast2smaller(rst,ren,cen,nonnan,sz)
% Which unsigned integer class  
[~,bin] = histc(numel(nonnan), [1 256 65536 4294967296 18446744073709551616]);
if bin == 5
    throwAsCaller(MException('MATLAB:pmaxsize','Maximum variable size allowed by the program is exceeded.'))
end

% Pick the conversion fun
fun = {@uint8, @uint16, @uint32, @uint64};
fun = fun{bin};

% Convert
rst    = fun(rst);
ren    = fun(ren);
cen    = fun(cen);
nonnan = fun(nonnan);
sz     = fun(sz);  
end