function s = pivotFromTo(tb)

% PIVOTFROMTO Pivots a table spreading the value between the from/to range
% 
%   ... = PIVOTFROMTO(TB) TB should be a table with:
%           - 1st col --> Id, i.e. the indicator variable
%           - 2nd col --> From date
%           - 3rd col --> To date
%           - 4rd col --> Value to pivot
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
s.Panel = unstack(tb, 'Val','Id');
fun = @(x) cast(cumsum(nan2zero(x)),classval);
s.Panel(:,2:end) = varfun(fun, s.Panel(:,2:end));

% Reset the one day extension of ending date
if yyyymmdd2serial(s.Panel.Date(end))-1 == yyyymmdd2serial(s.Panel.Date(end-1))
    s.Panel(end,:) = [];
else
    s.Panel.Date(end)  = serial2yyyymmdd(yyyymmdd2serial(s.Panel.Date(end))-1);
    s.Panel(end,2:end) = s.Panel(end-1,2:end);
end
s.Date = s.Panel.Date(2:end,1);
end