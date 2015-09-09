function s = pivotFromTo(tb, refdates)

% PIVOTFROMTO Pivots a table spreading the value between the from/to range
% 
%   PIVOTFROMTO(TB) TB should be a table with:
%           ID | From | To | Value
%       
%           The ID is the indexing var and Value is the unstacked variable.
%
%   PIVOTFROMTO(..., REFDATES) Retains the dates specified by REFDATES. 
%                              By default it keeps all days in the range
%                              [min(From), max(To)].
%
%   S = ... 
%       Is a structure with:
%           .Date   -->  unique vector of dates in ascending order
%           .Id     -->  unique vector of Ids in ascending order 
%           .Idname -->  name of the Id, e.g. 'permno'
%           .Panel  -->  unstacked table with pivoted values
%           
%   NOTE: all dates are uint32() in yyyymmdd format
%
% See also: UNSTACK, YYYYMMDD2SERIAL

oldnames = getVariableNames(tb);
tb       = setVariableNames(tb, {'Id','From','To','Val'});
classval = class(tb.Val);
tb.Val   = double(tb.Val);

% Extend To date by one day (will be nullified by the spread/pivoting)
to    = tb.To;
tb.To = serial2yyyymmdd(yyyymmdd2serial(to)+1);

% Stack start/end position to get +1/-1 Val
% NOTE: ensure that the original To is preserved
o  = zeros(size(to,1),1);
tb = table(repmat(tb.Id,3,1),...
           [tb.From; to;  tb.To],...
           [tb.Val;   o; -tb.Val],...
           'VariableNames',{'Id','Date','Val'});
tb = sortrows(tb,{'Date','Id'});

% Convert old id name capitalizing first letter only
s.Id     = unique(tb.Id);
s.Idname = lower(oldnames{1});

% Pivot and spread membership
s.Panel          = unstack(tb, 'Val','Id');
fun              = @(x) cumsum(nan2zero(x));
s.Panel(:,2:end) = tbextend.varfun(fun, s.Panel(:,2:end),'RenameVariables', false);

% Convert back to original class
if ~strcmpi(classval, 'double')
    s.Panel = convertColumn(s.Panel, classval, 2:size(s.Panel,2));
end

% Sample at reference dates
if nargin < 2
    s.Date = serial2yyyymmdd(yyyymmdd2serial(min(tb.Date)):yyyymmdd2serial(max(tb.Date)));
else
    s.Date = refdates;
end
s.Panel = sampledates(s.Panel,s.Date,true);
end