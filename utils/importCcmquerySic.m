function ccmsic = importCcmquerySic(path2zip)
% IMPORTCCMQUERYSIC Imports the zipped (CRSP) CCM webquery dataset with SIC codes into a table
%
%   IMPORTCCMQUERYSIC (PATH2ZIP)
%
% Missing codes (Ret):
%
%   CRSP code | SAS code | Description
%      -44         .E      Still valid
%
% See for details: https://wrds-web.wharton.upenn.edu/wrds/support/Data/_001Manuals%20and%20Overviews/_002CRSP/ccm-overview.cfm

zipfile = 'CRSPccmquery_sic.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
% csvfile = {'D:\TAQ\HFbetas\data\CRSP\55855a2eb0612cfd.csv'};

% Get header
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));
txt     = textscan(fid, '%u32%u16%*s%*s%*s%u32%*u32%u32%s','Delimiter',',');

% Convert 'E' end dates to 99999999
idx           = ismember(txt{5}, 'E');
enddate       = zeros(numel(idx),1,'uint32');
enddate(idx)  = 99999999;
tmp           = textscan(char(txt{5}(~idx))', '%8u32');
enddate(~idx) = tmp{1};
txt{5}        = enddate;

% Convert to table
ccmsic = table(txt{:}, 'VariableNames', headers([1:2,6,8:9]));

% Save
filename = sprintf('%s_ccmsic.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'ccmsic')

end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end