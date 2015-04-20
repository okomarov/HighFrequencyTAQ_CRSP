function tf = issp500member(tb)
% ISSP500MEMBER Checks which Id - Date pair are sp500 members
%
%   ISSP500MEMBER(TB) TB is a table with Id and yyyymmdd Date 

if isa(tb,'dataset')
    tb = dataset2table(tb);
end

mst = loadresults('issp500');
mst = mst(mst.Issp500,{'Permno','Date'});
tf  = ismembIdDate(tb.Permno, tb.Date, mst.Permno, mst.Date);
end