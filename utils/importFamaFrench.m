function FFfactors = importFamaFrench(path2zip)
% IMPORTFAMAFRENCH Imports the zipped Fama and French factors
%
%   IMOPRTMSENAMES(PATH2ZIP)
% 
% See <a href="http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html">Fama French data</a>

zipfile = 'F-F_Research_Data_Factors_daily.zip';
writeto = '.\results';

% Unzip and loadf as row strings
txtfile = unzip(fullfile(path2zip, zipfile),path2zip);
fid     = fopen(char(txtfile));
vnames  = textscan(fid, '%s%s%s%s',1,'MultipleDelimsAsOne',true,'HeaderLines',4);
txt     = textscan(fid, '%d%f%f%f%f','MultipleDelimsAsOne',true);
fclose(fid);

% Convert into table
FFfactors = table(txt{:},'VariableNames',['Date', 'MktMinusRF', vnames{2:end}]);

% Delete unzipped .txt
delete(txtfile{:})

% Save
fname = sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'FFfactors');
save(fullfile(writeto, fname), 'FFfactors')
end