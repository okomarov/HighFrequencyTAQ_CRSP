function res = importCrspIndices(fullzipfile)
% IMPORTCRSPINDICES Imports the zipped CRSP msix, dsix and similar datasets into a table
%
%   IMPORTCRSPINDICES(FULLZIPFILE) 
%       Full path to the zipfile

writeto = '.\results';

[~,datasetname] = fileparts(fullzipfile);

% Unzip and loadf as row strings
csvfile  = unzip(fullzipfile,tempdir());
fid      = fopen(char(csvfile));
txt      = textscan(fid, '%s',1,'Delimiter','');
txt      = txt{1};
final_cleanup = onCleanup(@()myCleanup(fid,csvfile{1})); 

% Headers
vnames = upperfirst(strsplit(txt{1},','));

% Parse fields
fmt  = ['%u32', repmat('%f',1,34)];
data = textscan(fid, fmt, 'HeaderLines',1,'Delimiter',',');

% Convert into table and replace missing
res = table(data{:},'VariableNames',vnames);

% Save
fname = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),datasetname);
save(fullfile(writeto,fname), 'res')
end

function myCleanup(fid,filename)
fclose(fid);
delete(filename);
end
