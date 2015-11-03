function tb = addPermno(tb)
% Add permno to table with Id and Date pairs
try
    tb.Permno(1);
    warning('addPermno:hasPermno','The table already has a PERMNO.')
    return
catch
end

try
    mst = loadresults('masterPermno');
catch
    mst = loadresults('masterPermno','..\results');
end
[idx,pos] = ismembIdDate(tb.Id,tb.Date,mst.Id, mst.Date);
if nnz(idx)
    warning('addPermno:matchFail','Found non-matching "Id - Date" pairs.')
end
tb.Permno(idx,1) = mst.Permno(pos(idx));
end