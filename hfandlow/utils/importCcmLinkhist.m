function ccm = importCcmLinkhist(path2zip)
% MIPORTCCMLINKHIST Imports the zipped CRSP/COMPUSTAT merged dataset into a table
%
%   MIPORTCCMLINKHIST(PATH2ZIP)
%
% Link types:
% 
%   LinkId | Linktype | Priority     | Description 
%      1        LC      Primary        Link research complete by CRSP. Standard connection between databases.
%      2        LU      Primary        Link is unresearched by CRSP. It is established by comparing the Compustat and historical CRSP CUSIPs. Most popular.
%      3        LS      Primary        Link valid for this security only. Other CRSP PERMNOs with the same PERMCO will link to other GVKEYs. LS links mainly relate to ETFs where a single CRSP PERMCO links to multiple Compustat GVKEYs
%      4        LX      Secondary      Link to a security that trades on foreign exchange not included in CRSP data.
%      5        LD      Secondary      Duplicate link to a security. Two GVKEYs map to a single PERMNO (PERMCO) during the same period, and this link should not be used. Almost all of these cases happened before 1990.
%      6        LN      Secondary      Primary link exists but Compustat does not have prices.
%      7        NP      Non-matching    
%      8        NR      Non-matching   No link available; confirmed by research.
%      9        NU      Non-matching   No link available; not yet confirmed.
%
% See for details: 
%   * https://wrds-web.wharton.upenn.edu/wrds/support/Data/_001Manuals%20and%20Overviews/_002CRSP/ccm-overview.cfm
%   * https://wrds-web.wharton.upenn.edu/wrds/E-Learning/_000Course%20Materials/Overview%20of%20CCM.pdf.cfm

zipfile = 'CRSPccmxpf_linkhist.csv.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile = unzip(fullfile(path2zip, zipfile),path2zip);
fid     = fopen(char(csvfile));
cleanup = onCleanup(@()cleanupFile(fid));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = upperfirst(regexp(headers{1}{1},',','split'));

% Parse fields
fmt  = '%u32%c%s%s%u32%u32%u32%s';
data = textscan(fid, fmt, 'Delimiter',',');

% Adjust end dates (E to 99999999)
idx           = strcmpi(data{end},'E');
enddate       = repmat(uint32(99999999), numel(idx), 1);
c             = [char(data{end}) repmat(' ',numel(idx),1)];
tmp           = textscan(c(~idx,:)', '%9u32');
enddate(~idx) = tmp{1};

% Convert to table
ccm           = table(data{:}, 'VariableNames', headers);
ccm.Linkenddt = enddate;

% Map Linktype
[~, ccm.Linktype] = ismember(ccm.Linktype, {'LC';'LU';'LS';'LX';'LD';'LN';'NP';'NR';'NU'});
ccm.Linktype      = uint8(ccm.Linktype);

% Save
filename = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'ccm');
save(fullfile(writeto, filename), 'ccm')
end

function cleanupFile(fid)
fname = fopen(fid);
fclose(fid);
delete(fname)
end