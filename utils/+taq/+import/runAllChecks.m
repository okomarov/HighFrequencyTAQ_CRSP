function runAllChecks(matfolder,testfolder)
% Day spillover
prob = checkDaySpillover(matfolder);
if ~isempty(prob)
    disp(prob)
    error('There are spillovers of daily data across multiple files.')
end

% Check master from-to continuity
dd         = dir(fullfile(matfolder,'*.idx'));
indexnames = {dd.name};
for ii = 1:numel(indexnames)
    load(fullfile(matfolder,indexnames{ii}),'-mat');
    if any(index.To(1:end-1)+1 ~= index.From(2:end))
        error('File %d has gaps between To and the next From.',ii)
    end
end
end