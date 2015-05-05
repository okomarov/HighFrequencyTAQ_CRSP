function s = pivotFromTo(tb)

% PIVOTFROMTO Pivots a table spreading the value between the from/to range
% 
%   PIVOTFROMTO(TB) TB should be a table with:
%           ID | From | To | Value
%       
%           The ID is the indexing var and Value is the unstacked variable.
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

% Reference dates
s.Date = union(tb.From, tb.To);

% Extend To date by one day (will be nullified by the spread/pivoting)
tb.To = serial2yyyymmdd(yyyymmdd2serial(tb.To)+1);

% Stack start/end position to get +1/-1 Val
tb = table(repmat(tb.Id,2,1),...
           [tb.From; tb.To],...
           [tb.Val; -tb.Val],...
           'VariableNames',{'Id','Date','Val'});
tb = sortrows(tb,{'Date','Id'});

% Convert old id name capitalizing first letter only
s.Id     = unique(tb.Id);
s.Idname = lower(oldnames{1});

% Pivot and spread membership
s.Panel          = tbextend.unstack(tb, 'Val','Id');
fun              = @(x) cumsum(nan2zero(x));
s.Panel(:,2:end) = tbextend.varfun(fun, s.Panel(:,2:end),'RenameVariables', false);
s.Panel          = convertColumn(s.Panel, classval, 2:size(s.Panel,2));

% Sample at reference dates
s.Panel = sampledates(s.Panel,s.Date,true);
end