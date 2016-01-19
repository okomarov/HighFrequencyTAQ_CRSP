function tb = addOvernightRet(tb, reton)
% Returns are treated as log returns

% See also: estimateOverIntraCloseRet

if nargin < 2
    try
        reton = loadresults('return_intraday_overnight');
    catch
        reton = loadresults('return_intraday_overnight','..\results');
    end
end

try
    tb.Date(1);
catch
    tb.Date = serial2yyyymmdd(tb.Datetime);
end

% Find where overnight should be positioned
iDateChange   = [true; diff(int32(tb.Date)) ~= 0];
iPermnoChange = [true; diff(int32(tb.Permno)) ~= 0];
posOn         = find(iPermnoChange | iDateChange);

% Check if has already a return
iadd      = ~isnan(tb.Ret(posOn));
[idx,pos] = ismembIdDate(tb.Permno(posOn), tb.Date(posOn),reton.Permno,reton.Date);

tb.Ret(posOn(idx &  iadd)) = tb.Ret(idx & iadd) + reton.RetCO(pos(idx & iadd));
tb.Ret(posOn(idx & ~iadd)) = reton.RetCO(pos(idx & ~iadd));
end