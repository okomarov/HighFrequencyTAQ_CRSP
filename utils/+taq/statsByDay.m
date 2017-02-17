function res = statsByDay(path2data, outname, fun)
% STATSBYDAY Wrapper for the +statsByDay package. Do NOT call directly!
res = iterate_datastore(path2data, outname, fun);
end
