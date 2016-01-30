function runAllChecks(matfolder,testfolder)
prob = checkDaySpillover(matfolder);
if ~isempty(prob)
    disp(prob)
    error('There are spillovers of daily data across multiple files.')
end
end