function names = importNames(path2zip)
% IMPORTNAMES Imports the zipped TAQ NAMES dataset into a table
%
%   IMPORTNAMES(PATH2ZIP) 

zipfile = 'NAMES.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile  = unzip(fullfile(path2zip, zipfile),path2zip);
fid      = fopen(char(csvfile));
cleanup  = onCleanup(@() myCleanupFunction(fid,csvfile));
txt      = textscan(fid, '%s','Delimiter','');
txt      = txt{1};

% Headers
vnames = setdiff(regexp(txt{1},'[^\w]*','split'),'','stable');

% Parse fields
txt  = strcat(txt,',');
fmt  = '%s%q%s%u32%u32%u16';
data = textscan([txt{2:end}], fmt, 'Delimiter',',');

% Convert into table and replace missing
names = table(data{:},'VariableNames',vnames);

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'TAQ_names')), 'names')
end
function myCleanupFunction(fid,csvfile)
fclose(fid);
delete(char(csvfile));
end