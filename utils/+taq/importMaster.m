function importMaster(rootfolder)
% IMPORTMST Imports .csv and .tab master files into one table
% 
%   IMPORTMST(PATH2ZIPS) The .csv and .tab files are both zipped in two
%                        archives that should reside on PATH2ZIPS
%
%   Format ot the output table
%
%       SYMBOL | NAME | CUSIP | ETN | ETA | ETB | ETP | ETX | ETT | ETO| ...
%       cstr   | cstr | cstr  | tf  | tf  | tf  | tf  | tf  | tf  | tf | ...
%
%       ... ETW | ITS | ICODE | SHROUT | UOT | DENOM | TYPE | FDATE
%       ... tf  | tf  | cstr  | double | u32 | char  |  i8  | u32

%% EDIT THESE
masterfolder = fullfile(rootfolder,'raw','master');
outdir       = fullfile(rootfolder,'master');
%% Import

% Manual CSV
% ==========
fprintf('Processing manual .csv files\n')
list   = unzip(fullfile(masterfolder,'TAQmast.csv.zip'),tempdir());
list   = sort(list);
% Due to a change in FDATE as per TAQ readme.txt which screwed the field,exclude 2000/02 and 2000/03 master files
idx    = cellfun('isempty', regexp(list,'(mst200002|mst200003)','once'));
delete(list{~idx})
list   = list(idx);
nfiles = numel(list);
cleanup = onCleanup(@() delete(list{:}));

for ii = 1:nfiles
    disp(ii/nfiles*100)
    taq.import.masterFile(list{ii}, outdir)
end

% Manual TAB
% ==========
fprintf('Processing .tab files\n')
list   = unzip(fullfile(masterfolder,'TAQmast.tab.zip'),tempdir());
list   = sort(list);

% Only between 2010/07 and 2010/12
idx    = ~cellfun('isempty', regexp(list,'2010','once'));
delete(list{~idx})
list   = list(idx);
nfiles = numel(list);
cleanup = onCleanup(@() delete(list{:}));

opt.UseTextscan = true;
opt.ImportFmt   = '%10c %30c %12c %c%c%c %*2c %c%c%c%c%c %c %4c %10c %4c %c %c %8c %*[^\n]';
opt.VarNames    = {'SYMBOL','NAME','CUSIP','ETN','ETA','ETB','ETP','ETX','ETT','ETO','ETW','ITS','ICODE','SHROUT','UOT','DENOM','TYPE','FDATE'};
opt.ImportOther = {'Whitespace',''};

for ii = 1:nfiles
    disp(ii/nfiles*100)
    taq.import.masterFile(list{ii}, outdir, opt)
end

% Automated CSV
% =============
fprintf('Processing automated .csv files\n')
ziplist = dir(fullfile(masterfolder,'MAST_*.zip'));
ziplist = sort(fullfile(masterfolder,{ziplist.name}))';
nfiles  = numel(ziplist);

for ii = 1:nfiles
    disp(ii/nfiles*100)
    fname = unzip(ziplist{ii},tempdir());
    cleanup = onCleanup(@() delete(fname{1}));
    taq.import.masterFile(fname{1}, outdir)
end
end
