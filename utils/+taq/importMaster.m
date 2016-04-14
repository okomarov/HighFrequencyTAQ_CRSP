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
fprintf('Unzipping manual .csv files\n')
list   = unzip(fullfile(masterfolder,'TAQmast.csv.zip'),tempdir());
list   = sort(list);
% Due to a change in FDATE as per TAQ readme.txt which screwed the field,exclude 2000/02 and 2000/03 master files
idx    = cellfun('isempty', regexp(list,'(mst200002|mst200003)','once'));
delete(list{~idx})
list   = list(idx);
nfiles = numel(list);

for ii = 1:nfiles
    taq.import.masterFile(list{ii}, outdir)
end
delete(list{:})

% Manual TAB
% ==========
fprintf('Unzipping .tab files\n')
list   = unzip(fullfile(masterfolder,'TAQmast.tab.zip'),tempdir());
list   = sort(list);

% Only between 2010/07 and 2010/12
idx    = ~cellfun('isempty', regexp(list,'2010','once'));
delete(list{~idx})
list   = list(idx);
nfiles = numel(list);

opt.UseTextscan = true;
opt.ImportFmt   = '%10c %30c %12c %c%c%c %*2c %c%c%c%c%c %c %4c %10c %4c %c %c %8c %*[^\n]';
opt.VarNames    = {'SYMBOL','NAME','CUSIP','ETN','ETA','ETB','ETP','ETX','ETT','ETO','ETW','ITS','ICODE','SHROUT','UOT','DENOM','TYPE','FDATE'};
opt.ImportOther = {'Whitespace',''};

for ii = 1:nfiles
    disp(ii)
    taq.import.masterFile(list{ii}, outdir, opt)
end
delete(list{:})

% Automated CSV
% =============
ziplist = dir(fullfile(masterfolder,'MAST_*.zip'));
ziplist = sort(fullfile(masterfolder,{ziplist.name}))';
nfiles  = numel(ziplist);

for ii = 1:nfiles
    disp(ii)
    fname = unzip(ziplist{ii},tempdir());
    taq.import.masterFile(fname{1}, outdir)
    delete(fname{1})
end
end
