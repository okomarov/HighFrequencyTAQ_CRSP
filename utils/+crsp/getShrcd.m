function shrcd = getShrcd()
% GETSHRCD Loads share code for Permno and Date pairs

folder = fullfile(fileparts(mfilename('fullpath')),'data');
try
    shrcd = loadlatest('shrcd',folder);
catch
    warning('Making share code dataset.')
    shrcd = makeShrcd(folder);
end
end

function shrcd = makeShrcd(folder)
try
    msenames = loadlatest('msenames',folder);
catch
    error('No MSENAMES dataset found. Import that first.')
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
msenames = stack(msenames.Panel, vnames(2:end),...
                'IndexVariableName', 'Permno',...
                'NewDataVariableName','Shrcd');

% Keep only non null values
msenames = msenames(msenames.Shrcd ~= 0,:);

% Convert xPermno to numeric permno
msenames.Permno = uint32(xstr2num(msenames.Permno));

% Save
shrcd = msenames;
fname = fullfile(folder,sprintf('shrcd_%s.mat',datestr(now,'yyyymmdd_HHMM')));
save(fname, 'shrcd')
end
