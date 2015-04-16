function shrcd = mapShrcd2mst

% MAPSHRCD2MST Maps share code type to master records with Id and Date and Symbol

% Load master permno
mst = loadresults('masterPermno');

% Load msenames
try
    msenames = loadresults('msenames');
catch
    msenames = importMsenames('.\data\CRSP\');
end
msenames             = msenames(msenames.NAMEENDT > 19923112,:);
idx                  = msenames.NAMEDT <= 19923112;
msenames.NAMEDT(idx) = 19930101;

% Filter by permno
idx      = ismember(msenames.PERMNO, unique(mst.Permno));
msenames = msenames(idx,:);

% Time consolidation
msenames = sortrows(msenames(:,{'PERMNO','NAMEDT','SHRCD'}));
idx      = isfeatchange(msenames(:,{'PERMNO','SHRCD','NAMEDT'}));
msenames = msenames(idx,:);

% Unstack
msenames = convertColumn(msenames,'double','SHRCD');
msenames = unstack(msenames,'SHRCD','PERMNO');

% Sample dates
msenames = renameVarNames(msenames, 'Date','NAMEDT');
msenames = sampledates(msenames, unique(mst.Date));

% Stack back into tall table
vnames   = getVariableNames(msenames);
msenames = stack(msenames, vnames(2:end), ...
                'IndexVariableName', 'Permno',...
                'NewDataVariableName','Shrcd');

% Keep only non 0 values
msenames = msenames(msenames.Shrcd ~= 0,:);

% Convert literal xPermno to numeric id
Permno          = char(msenames.Permno);
Permno          = Permno(:,2:end)';
Permno          = textscan([Permno; repmat(' ',1,size(Permno,2))],'%u32');
msenames.Permno = Permno{1};

% Map to mst
[idx,pos]        = ismembIdDate(mst.Permno,mst.Date, msenames.Permno, msenames.Date);
mst.Shrcd(idx,1) = uint8(msenames.Shrcd(pos(idx)));

% Save
shrcd = mst;
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'shrcd')), 'shrcd')

end