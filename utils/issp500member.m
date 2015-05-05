function tf = issp500member(tb)
% ISSP500MEMBER Checks which Id - Date pair are sp500 members
%
%   ISSP500MEMBER(TB) TB is a table with Id and yyyymmdd Date 

if isa(tb,'dataset')
    tb = dataset2table(tb);
end

spconst = loadresults('spconst');
tf      = ismembIdDate(tb.Permno, tb.Date, spconst.Permno, spconst.Date);
end