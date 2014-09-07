function shrcd = mapShrcd2mst

% MAPSHRCD2MST Maps share code type to master records with Id and Date and Symbol

% Load msenames
try
    loadresults('msenames')
catch
    msenames = importMsenames('.\data\CRSP\');
end
% Sort
msenames = sortrows(msenames(:,{'PERMNO','NAMEDT','SHRCD'}));
% Time consolidation
idx      = isfeatchange(msenames(:,{'PERMNO','SHRCD','NAMEDT'}));
msenames = msenames(idx,:);

% Load taq2crsp
loadresults('taq2crsp')
if isa(taq2crsp,'dataset'), taq2crsp = dataset2table(taq2crsp); end
taq2crsp = taq2crsp(~isnan(taq2crsp.permno),{'ID','permno','datef'});
% Time consolidation
taq2crsp = sortrows(taq2crsp,{'ID','datef'});
idx      = isfeatchange(taq2crsp);
taq2crsp = taq2crsp(idx,:);

% Unstack msenames into shrcd
msenames       = msenames(ismember(msenames.PERMNO, taq2crsp.permno),:);
msenames.SHRCD = single(msenames.SHRCD); % For NaN padding in unstack
shrcd          = unstack(msenames,'SHRCD','PERMNO');
shrcd          = sortrows(shrcd,'NAMEDT');
% Fill gaps
names          = getVariableNames(shrcd);
shrcd(:,2:end) = varfun(@nanfillts,shrcd(:,2:end));
shrcd          = setVariableNames(shrcd,['Date', names(2:end)]);

% Unstack taq2crsp into mask
id2permno = unique(taq2crsp(:,{'ID','permno'}));
id2permno.ID = uint16(id2permno.ID);
import matlab.lang.*
id2permno.Properties.RowNames = makeValidName(cellstr(num2str(id2permno.ID)));
taq2crsp.Val = ones(size(taq2crsp,1),1,'single');
mask     = unstack(taq2crsp(:,{'Val','ID','datef'}),'Val','ID');
mask     = sortrows(mask,'datef');
mask.Properties.VariableNames{1} = 'Date';

% Sample dates
loadresults('uniqueID')
refdates = unique(uniqueID.Date);
shrcd    = sampledates(shrcd,refdates);
mask     = sampledates(mask,refdates);

% Create mask 
mask = table2array(mask(:,2:end));
mask = logical(uint8(nanfillts(mask)));

% Sample shrcd from msenames
permnos2get = makeValidName(cellstr(num2str(id2permno.permno)));
shrcd = uint8(table2array(shrcd(:, permnos2get)));

% Apply mask 
shrcd(~mask) = 0; clear mask

shrcd = array2table(shrcd,'VariableNames',id2permno.Properties.RowNames);
shrcd.Date = refdates;

% Stack back
shrcd = stack(shrcd,id2permno.Properties.RowNames,'NewDataVariableName','Shrcd',...
              'IndexVariableName','UnID');

% Clean up
shrcd = shrcd(shrcd.Shrcd ~= 0,:);
[~,pos] = ismember(shrcd.UnID,id2permno.Properties.RowNames);
shrcd.UnID = id2permno.ID(pos); clear pos

% Map to uniqueID (in blocks)
uniqueID.Shrcd =  zeros(size(uniqueID,1),1,'uint8');
N = 4;
edges = fix(linspace(0,size(shrcd,1),N+1));
for ii = 1:N
    [idx,pos] = ismember(uniqueID(:,{'Date','UnID'}),shrcd(edges(ii)+1:edges(ii+1),{'Date','UnID'}));
    uniqueID.Shrcd(idx) = shrcd.Shrcd(pos(idx)+edges(ii));
end

% Save
shrcd = uniqueID;
save(fullfile('.\results',sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'shrcd')), 'shrcd')

end