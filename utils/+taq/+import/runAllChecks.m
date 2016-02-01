function runAllChecks(matfolder,testfolder)
% Day spillover
prob = checkDaySpillover(matfolder);
if ~isempty(prob)
    disp(prob)
    error('There are spillovers of daily data across multiple files.')
end

% Check master from-to continuity
dd       = dir(fullfile(matfolder,'*.mst'));
mstnames = {dd.name};
for ii = 1:numel(mstnames)
    load(fullfile(path2data,mstnames{ii}),'-mat');
    if any(mst.From(2:end)-mst.To(1:end-1) ~= 1)
        error('File %d has gaps between To and the next From.',ii)
    end
end
end