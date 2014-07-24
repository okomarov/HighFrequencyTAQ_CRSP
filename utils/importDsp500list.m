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

% Add value
spconst.Val = true(size(spconst,1),1);

% Pivot
spconst = pivotFromTo(spconst);

% Save
save(fullfile(writeto,sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'spconst')), 'spconst')
end