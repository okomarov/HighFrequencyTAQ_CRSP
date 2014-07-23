function dseshares = importDseshares(path2zip)
% IMPORTDSESHARES Imports the zipped CRSP dseshares dataset into a table
%
%   IMPORTDSESHARES(PATH2ZIP) 

zipfile = 'CRSPdseshares.csv.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile  = unzip(fullfile(path2zip, zipfile),path2zip);
fid      = fopen(char(csvfile));
txt      = textscan(fid, '%s','Delimiter','');
txt      = txt{1};
fclose(fid);

% Headers
vnames = setdiff(regexp(txt{1},'[^A-Z_]*','split'),'','stable');

% Parse fields
fmt  = '%u32 %u32 %u32 %u32 %u8 %u16 %u32 %u32 %u8 %u16 %8c';
data = textscan([txt{2:end}], fmt, 'Delimiter',',');

% Convert into table and replace missing
dseshares = table(data{:},'VariableNames',vnames);

% Delete unzipped .csv
delete(csvfile{:})

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'dseshares')), 'dseshares')
end