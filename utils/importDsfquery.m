function dsfquery = importDsfquery(path2zip)
% IMPORTDSFQUERY Imports the zipped CRSP dsfquery dataset into a table
%
%   IMPORTDSFQUERY (PATH2ZIP)
%
% Missing codes (Ret):
%
%   CRSP code | SAS code | Description
%      -44         .E      No valid comparison for an excess return
%      -55         .D      No listing information
%      -66         .C      No valid previous price or > 10 periods before
%      -77         .B      Not trading on the current exchange
%      -88         .A      No data available to calculate returns
%      -99         .       No valid price (usually suspension or trading on unknown exchange)
%
% See for details: https://wrds-web.wharton.upenn.edu/wrds//query_forms/variable_descriptions.cfm?legacyPath=/ds/crsp/stock_a&legacyFile=msf.cfm&title=CRSP%20Annual%20Stock&libraryCode=crspa&fileCode=msf&adminEditPath=/wrds/ds/crsp/stock_a/msf_doc.cfm

zipfile = 'CRSPdsfquery.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
% csvfile = {'D:\TAQ\HFbetas\data\CRSP\3a010f6b509a2610.csv'};

% Get header
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

% Block process (ret needs additional processing)
dsfquery     = [];
missingcodes = [];
c            = 0;
N            = 1e6;
while ~feof(fid)
    % Parse c-th block
    c = c+1;
    disp(c)
    
    % Note that Returns might have char missing code
    txt = textscan(fid, '%u32%u32%f32%u64%s%f32%f32',N,'Delimiter',',');
    
    % Import delisting returns and returns
    [txt{5}, missing{2}] = dealWithMissingCodes(txt{5});
    
    % Expand pre-allocation
    if mod(c,100) == 1
        dsfquery = [dsfquery; cell(100,1)];
    end
    % Convert into table
    dsfquery{c} = table(txt{:}, 'VariableNames',headers(:));
end

% Concatenate
dsfquery = cat(1,dsfquery{:});

% Adjust for delisting returns as in Beaver, McNichols, Price 2007
dsfquery = adjustDelist(dsfquery);

% Save
filename = sprintf('%s_dsfquery.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'dsfquery')
end

function [out, codes] = dealWithMissingCodes(c)
% Find empty returns
iempty = cellfun('isempty',c);

% Find missing codes in returns and map them out
imisscode            = false(size(iempty));
imisscode(~iempty)   = ~cellfun('isempty',regexp(c(~iempty), '[ABCDEPST]','once'));
codes                = uint8(imisscode);
[~,codes(imisscode)] = ismember(c(imisscode), {'A';'B';'C';'D';'E';'P';'S';'T'});

% Convert into numbers
out      = NaN(size(c));
idx      = ~(iempty | imisscode);
tmp      = textscan(char(c(idx))', '%9f64');
out(idx) = tmp{1};
end

function dsf = adjustDelist(dsf)
% Load delisting returns
try
    delist = loadresults('dsedelist');
catch
    delist = importDsedelist('dsedelist');
end
% Estimate replacement values by delisting code for delisting dates with 
% missing delisting returns
[unid, ~, subs]    = unique(delist.Dlstcd);
repval             = accumarray(subs, delist.Dlret,[],@nanmean);
inan               = isnan(delist.Dlret);
delist.Dlret(inan) = repval(subs(inan));

% Adjust returns in dsfquery
[idx,pos]    = ismembIdDate(dsf.Permno,dsf.Date,delist.Permno,delist.Dlstdt);
dsf.Ret(idx) = (1 + dsf.Ret(idx)) .* (1 + delist.Dlret(pos(idx))) - 1;
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end