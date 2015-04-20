function tf = iscommonshare(tb)
% ISCOMMONSHARE Checks which Id - Date pairs are common shares (share type code 10 and 11) 
%
%   ISCOMMONSHARE(TB) TB is a table with Id and yyyymmdd Date 
if isa(tb,'dataset')
    tb = dataset2table(tb);
end

shrcd = loadresults('shrcd');
shrcd = shrcd(shrcd.Shrcd == 11 | shrcd.Shrcd == 10,{'Permno','Date'});
tf    = ismembIdDate(tb.Permno, tb.Date, shrcd.Permno, shrcd.Date);

end