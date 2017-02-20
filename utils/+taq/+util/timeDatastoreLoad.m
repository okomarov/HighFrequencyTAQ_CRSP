function res = timeDatastoreLoad(path2data)
res = taq.iterate_datastore(path2data, [], @(s,cached,filenum) void_());
end

function out = void_()
out = [];
end
