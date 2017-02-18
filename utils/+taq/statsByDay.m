function res = statsByDay(path2data, outname, fun, iterOpts)
% STATSBYDAY Wrapper for the +statsByDay package. Do NOT call directly!
iterOpts = [fieldnames(iterOpts), struct2cell(iterOpts)]';
res      = iterate_datastore(path2data, outname, fun, iterOpts{:});
end
