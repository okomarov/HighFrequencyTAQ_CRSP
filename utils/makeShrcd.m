function shrcd = makeShrcd()
% Load msenames
try
    msenames = loadresults('msenames');
catch
    msenames = importMsenames('.\data\CRSP\');
end
msenames             = msenames(:,{'PERMNO','NAMEDT','NAMEENDT','SHRCD'});
msenames             = msenames(msenames.NAMEENDT >= 19930101,:);
idx                  = msenames.NAMEDT < 19930101;
msenames.NAMEDT(idx) = 19930101;

% Time consolidation
msenames          = sortrows(msenames,{'PERMNO','NAMEDT'});
idx               = isfeatchange(msenames(:,{'PERMNO','SHRCD','NAMEDT'}));
st                = find(idx);
en                = [st(2:end)-1; size(msenames,1)];
enddate           = msenames.NAMEENDT(en);
msenames          = msenames(st,:);
msenames.NAMEENDT = enddate;

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