function shrcd = makeShrcd()
% Load msenames
try
    msenames = loadresults('msenames');
catch
    msenames = importMsenames('.\data\CRSP\');
end
msenames             = msenames(:,{'Permno','Namedt','Nameendt','Shrcd'});
msenames             = msenames(msenames.Nameendt >= 19930101,:);
idx                  = msenames.Namedt < 19930101;
msenames.Namedt(idx) = 19930101;

% Time consolidation
msenames = sortrows(msenames,{'Permno','Namedt'});
msenames = consolidateFromTo(msenames);

% Pivot ranges 
msenames = pivotFromTo(msenames);

% Stack back to tall
vnames   = getVariableNames(msenames.Panel);
msenames = stack(msenames.Panel, vnames(2:end), ...
                'IndexVariableName', 'Permno',...
                'NewDataVariableName','Shrcd');

% Keep only non null values
msenames = msenames(msenames.Shrcd ~= 0,:);

% Convert xPermno to numeric permno
msenames.Permno = uint32(xstr2num(msenames.Permno));

% Save
shrcd = msenames;
fname = fullfile('.\results',sprintf('%s_shrcd.mat',datestr(now,'yyyymmdd_HHMM')));
save(fname, 'shrcd')
end