function msenames = importMsenames(path2zip)
% IMPORTMSENAMES Imports the zipped CRSP msenames dataset into a table
%
%   IMOPRTMSENAMES(PATH2ZIP)

zipfile = 'CRSPmsenames.csv.zip';
writeto = '.\results';

% Unzip and loadf as row strings
csvfile  = unzip(fullfile(path2zip, zipfile),path2zip);
fid      = fopen(char(csvfile));
txt      = textscan(fid, '%s','Delimiter','');
fclose(fid);

% Replace empty and last column
txt = regexprep(txt{:},'""','"N/A"');
txt = regexprep(txt,',0$',char(10));

% Headers
vnames = setdiff(regexp(txt{1},'[^A-Z_]*','split'),'','stable');

% Parse fields
fmt    = ['"%u32" "%d32" "%d32" "%u8" "%d8" "%u16" %s %s %s %s '...
          '%s "%u32" %s %s %s "%u32" "%u32" "%u32" "%u8" "%u16" %s%*[\n]'];
data   = textscan([txt{2:end}], fmt, 'Delimiter',',','TreatAsEmpty','N/A');

% Convert dates to yyyymmdd fmt
c        = datenum('19600101','yyyymmdd');
f        = @(x) uint32(serial2yyyymmdd(double(x + c)));
data{2}  = f(data{2});
data{3}  = f(data{3});

% Get rid of the ""
icstr       = cellfun(@iscellstr, data);
data(icstr) = cellfun(@(x) regexprep(x,'"',''),data(icstr),'un',0);

% Convert into table and replace missing
msenames = table(data{:},'VariableNames',vnames(1:end-1));
msenames = standardizeMissing(msenames,'N/A');

% Delete unzipped .csv
delete(csvfile{:})

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'msenames')), 'msenames')
end