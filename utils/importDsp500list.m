function spconst = importDsp500list(path2csv)
% IMPORTDSP500LIST Imports the .csv CRSP dsp500list dataset into a structure
%
%   IMPORTDSP500LIST(PATH2CSV)

csvfile = 'dsp500list.csv';
writeto = '.\results';

% Load
spconst = readtable(fullfile(path2csv,csvfile),'Delimiter',',','ReadVariableNames',1);

% Pivot ranges 
spconst.Val = ones(size(spconst,1),1);
spconst     = pivotFromTo(spconst);

% Stack back to tall
vnames  = getVariableNames(spconst.Panel);
spconst = stack(spconst.Panel, vnames(2:end), ...
                'IndexVariableName', 'Permno',...
                'NewDataVariableName','Issp500');

% Keep only SP500 members
spconst         = spconst(spconst.Issp500 ~= 0,:);
spconst.Issp500 = [];

% Convert xPermno to numeric permno
spconst.Permno = uint32(xstr2num(spconst.Permno));

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'spconst')), 'spconst')
end