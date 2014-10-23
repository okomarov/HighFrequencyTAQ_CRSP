function dsfquery = importDsfquery(path2zip)
% IMPORTDSFQUERY Imports the zipped CRSP dsfquery dataset into a table
%
%   IMPORTDSFQUERY (PATH2ZIP) 
%
% Missing codes:
%
%   CRSP code | SAS code | Description
%      -44         .E      No valid comparison for an excess return
%      -55         .D      No listing information
%      -66         .C      No valid previous price or > 10 periods before
%      -77         .B      Not trading on the current exchange
%      -88         .A      No data available to calculate returns
%      -99         .       No valid price (usually suspension or trading on unknown exchange)

% See for details: https://wrds-web.wharton.upenn.edu/wrds//query_forms/variable_descriptions.cfm?legacyPath=/ds/crsp/stock_a&legacyFile=msf.cfm&title=CRSP%20Annual%20Stock&libraryCode=crspa&fileCode=msf&adminEditPath=/wrds/ds/crsp/stock_a/msf_doc.cfm

zipfile = 'CRSPdsfquery.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile  = unzip(fullfile(path2zip, zipfile),path2zip);

% Get header
fid     = fopen(char(csvfile));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = regexp(headers{1}{1},',','split');
headers = [upperfirst(headers), 'MissingId'];

% Block process (ret needs additional processing)
dsfquery = [];
c        = 0;
N        = 1e6;
while ~feof(fid)
    % Parse c-th block
    c  = c+1;
    disp(c)
    txt = textscan(fid, '%u32%u32%f32%f32%f32%s%f32',N,'Delimiter',',');
    
    % Find empty
    iempty    = cellfun('isempty',txt{6});
    
    % Find/map missing code
    imisscode = ~cellfun('isempty',regexp(txt{6}, '[ABCDE]','once'));
    code      = uint8(imisscode);
    [~,code(imisscode)]  = ismember(txt{6}(imisscode), {'A';'B';'C';'D';'E'});    
    
    % Fill them with 0s, parse into number and replace with NaNs
    idx         = iempty | imisscode;
    txt{6}(idx) = {'0'};
    txt(6)      = textscan(char(txt{6})', '%9f64');
    txt{6}(idx) = NaN;
    
    % Expand pre-allocation
    if mod(c,100) == 1
        dsfquery = [dsfquery; cell(100,1)];
    end
    % Convert into table
    dsfquery{c} = table(txt{:}, code, 'VariableNames',headers);
end
fclose(fid);

% Concatenate
dsfquery = cat(1,dsfquery{:});

% Delete unzipped .csv
delete(csvfile{:})

% Save
filename = sprintf('%s_dsfquery.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'dsfquery')

end