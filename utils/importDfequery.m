function dfequery = importDfequery(path2zip)
% IMPORTDFEQUERY Imports the zipped CRSP dfequery dataset into a table
%
%   IMPORTDFEQUERY (PATH2ZIP) 

%% Subs missing return codes and 'B','C' with NaN but add additional uint8 column to map it
%{
-66.0 Valid current price but no valid previous price. Either first price, unknown
exchange between current and previous price, or more than 10 periods
between time t and the time of the preceding price t’.
-77.0 Not trading on the current exchange at time t.
-88.0 No data available to calculate returns.
-99.0 Missing return due to missing price at time t; usually due to suspension
in trading or trading on unknown exchange.
%}

zipfile = 'CRSPdfequery.csv.zip';
writeto = '.\results';

% Unzip and load as row strings
csvfile  = unzip(fullfile(path2zip, zipfile),path2zip);

% Get header
fid     = fopen(char(csvfile));
headers = textscan(fid, '%s',1,'Delimiter','');
headers = regexp(headers{1}{1},',','split');

% Block process (ret needs additional processing)
dfequery = [];
c        = 0;
N        = 1e6;
while ~feof(fid)
    % Parse c-th block
    c  = c+1;
    disp(c)
    txt = textscan(fid, '%u32%u32%f32%f32%f32%s%f32',N,'Delimiter',',');
    
    % Find empty or unwanted chars in return field
    idx = ~cellfun('isempty',regexp(txt{1,6}, '[^0-9\.-]','once')) | cellfun('isempty',txt{1,6});
%     unique(txt{1,6}(idx,:))
    
    % Fill them with 0s, parse into number and replace with NaNs
    txt{1,6}(idx) = {'0'};
    txt(1,6)      = textscan(char(txt{1,6})', '%9f64');
    
    % Expand pre-allocation
    if mod(c,100) == 1
        dfequery = [dfequery; cell(100,1)];
    end
    % Convert into table
    dfequery{c} = table(txt{:},'VariableNames',headers);
end
fclose(fid);

% Concatenate
dfequery = cat(1,dfequery{:});

% Delete unzipped .csv
delete(csvfile{:})

% Save
filename = sprintf('%s_dfequery.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile(writeto, filename), 'dfequery')

end