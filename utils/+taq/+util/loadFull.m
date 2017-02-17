function s = loadFull(matname)
load(matname);
load(strrep(matname,'.mat','.idx'),'-mat');
s = struct('data',{data},'index',{index},'symbol',{symbol});
end
