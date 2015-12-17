function splitMatTrades(path2data,nrecords)
% Read .mat filenames
mat = dir(fullfile(path2data,'*.mat'));
mst = dir(fullfile(path2data,'*.mst'));

% Preallocate
nfiles = numel(mat);
for f = 1:nfiles
    disp(f)
    sd = load(fullfile(path2mst,mat(f).name));
    sm = load(fullfile(path2mst,mst(f).name),'-mat');
    [~,~,bin] = histcounts(sm.mst.To,0:nrecords:size(sd.data,1));
    
end
    

end