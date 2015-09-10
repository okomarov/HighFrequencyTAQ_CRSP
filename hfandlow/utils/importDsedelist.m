function dsedelist = importDsedelist(path2zip)
% IMPORTDSEDELIST Imports the zipped CRSP dsedelist dataset into a table
%
%   IMPORTDSEDELIST (PATH2ZIP)
%
% Missing codes (Ret):
%
%   CRSP code | SAS code | Description
%     -55       .S         CRSP has no source to establish a value after delisting
%     -66       .T         More than 10 trading periods between a security's last price and its first available price on a new exchange
%     -88       .A         Security is still active
%     -99       .P         Security trades on a new exchange after delisting, but CRSP currently has no sources to gather price information
%
% See for details: https://wrds-web.wharton.upenn.edu/wrds//query_forms/variable_descriptions.cfm?legacyPath=/ds/crsp/stock_a&legacyFile=msf.cfm&title=CRSP%20Annual%20Stock&libraryCode=crspa&fileCode=msf&adminEditPath=/wrds/ds/crsp/stock_a/msf_doc.cfm

zipfile = 'CRSPdsedelist.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
% csvfile = {'D:\TAQ\HFbetas\data\CRSP\CRSPA.DSEDELIST.csv'};

% Get header
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));
txt     = textscan(fid, '%u32%u32%u16%*s%*s%*s%*s%*s%*s%*s%s%*[^\n]','Delimiter',',');

% Deal with missing codes
[txt{4}, missing{1}] = dealWithMissingCodes(txt{4});

% Convert to table
dsedelist = table(txt{:}, 'VariableNames', headers([1:3,11]));

% Save
filename = sprintf('%s_dsedelist.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'dsedelist')
end

function [out, codes] = dealWithMissingCodes(c)
% Find empty returns
iempty = cellfun('isempty',c);

% Find missing codes in returns and map them out
imisscode            = false(size(iempty));
imisscode(~iempty)   = ~cellfun('isempty',regexp(c(~iempty), '[APST]','once'));
codes                = uint8(imisscode);
[~,codes(imisscode)] = ismember(c(imisscode), {'A';'P';'S';'T'});

% Convert into numbers
out      = NaN(size(c));
idx      = ~(iempty | imisscode);
tmp      = textscan(char(c(idx))', '%9f64');
out(idx) = tmp{1};
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end