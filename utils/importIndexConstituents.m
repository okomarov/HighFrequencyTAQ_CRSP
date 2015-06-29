function [constit, desc] = importIndexConstituents(path2zip)
% IMPORTINDEXCONSTITUENTS Imports the zipped COMPUSTAT index constituents web-query into a table
%
%   IMPORTINDEXCONSTITUENTS(PATH2ZIP)

zipfile = 'COMPUSTATindexConstQuery.csv.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

% Parse fields
fmt  = '%u32%u32%u32%u32%s%s%s%s%s';
data = textscan(fid, fmt, 'Delimiter',',');

% Convert to table
constit = table(data{:}, 'VariableNames', headers);

% Factor out index info
desc    = unique(constit(:,[2,5:8]));
constit = constit(:,[1:4,9]);

% Adjust Thru dates (0 to 99999999)
idx               = constit.Thru == 0;
constit.Thru(idx) = 99999999;

% Save
filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'constituents');
save(fullfile(writeto, filename), 'constit','desc')
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end