function dsedist = importDsedist(path2zip)
% IMPORTDSEDIST Imports the zipped CRSP dsedist dataset into a table
%
%   IMPORTDSEDIST (PATH2ZIP)
%
%
% Note: total return is calculated by re-investing dividends of the
%       ex-dividend date (EXDT) according to the following formula:
%
%       r_t = (P_t*f_t + div_t)/P_t-1 - 1
%       
%       See CRSP's manual under "CRSP calculations > Return" on page 98.

zipfile = 'CRSPdsedist.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
% csvfile = {'D:\TAQ\HFbetas\data\CRSP\CRSPA.DSEDELIST.csv'};

% Get header
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

% Import
txt = textscan(fid, '%u32%u16%f%f%*f%u32%u32%*u32%u32%*[^\n]','Delimiter',',');

% Convert to table
dsedist = table(txt{:}, 'VariableNames', headers([1:4,6,7,9]));

% Save
filename = sprintf('%s_dsedist.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'dsedist')
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end