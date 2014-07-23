function spconst = importDsp500list(path2csv)
% IMPORTDSP500LIST Imports the .csv CRSP dsp500list dataset into a structure
%
%   IMPORTDSP500LIST(PATH2CSV)

aliveafter = 19921231;
csvfile    = 'dsp500list.csv';
writeto    = '.\results';

% Load
spconst = readtable(fullfile(path2csv,csvfile),'Delimiter',',','ReadVariableNames',1);
spconst = spconst(spconst.ending > aliveafter,:);

% Extend ending by one day (will be nullified by the spread/pivoting)
spconst.ending = serial2yyyymmdd(yyyymmdd2serial(double(spconst.ending))+1);

% Stack start/end position to get +1/-1 Val
nrows   = size(spconst,1);
spconst = table(repmat(spconst.PERMNO,2,1),...
    [spconst.start; spconst.ending],...
    [ones(nrows,1); -ones(nrows,1)],...
    'VariableNames',{'Permno','Date','Val'});
spconst = sortrows(spconst,{'Date','Permno'});

% Pivot and spread membership
panel   = unstack(spconst, 'Val','Permno');
dates   = panel.Date;

permnos = str2double(regexprep(panel.Properties.VariableNames(2:end),'x',''));
panel   = table2array(panel(:,2:end));
panel(isnan(panel)) = 0;
panel   = logical(cumsum(panel));

% Reset the one day extension of ending date
dates(end)   = serial2yyyymmdd(yyyymmdd2serial(double(dates(end)))-1);
panel(end,:) =  true;

% Collect
spconst = struct('Panel',panel,'Date',dates,'Permno',permnos);

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'spconst')), 'spconst')
end