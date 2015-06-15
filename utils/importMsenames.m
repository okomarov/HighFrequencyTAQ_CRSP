function msenames = importMsenames(path2zip)
% IMPORTMSENAMES Imports the zipped CRSP msenames dataset into a table
%
%   IMOPRTMSENAMES(PATH2ZIP)

zipfile = 'CRSPmsenames.csv.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

% Parse fields
fmt  = ['%u32%u32%u32%u8%u8%u16%s%s%s%s '...
          '%s%u32%s%s%s%u32%u32%u32%u8%u16%s'];
data = textscan(fid, fmt, 'Delimiter',',');

% Convert to table
msenames = table(data{:}, 'VariableNames', headers);

% Save
filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'msenames');
save(fullfile(writeto, filename), 'msenames')
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end