function TAQmaster = importTaqMaster(path2zips)
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
csvtables = 'TAQmast.csv.zip';
tabtables = 'TAQmast.tab.zip';
writeto   = '.\results';

%% Load .csv
% List all TAQ .csv master files
fprintf('Unzipping .csv files\n')
list      = unzip(fullfile(path2zips,csvtables),path2zips);
list      = sort(list);
% Due to a change in FDATE as per TAQ readme.txt which screwed the field,exclude 2000/02 and 2000/03 master files
idx       = cellfun('isempty', regexp(list,'(mst200002|mst200003)','once'));
delete(list{~idx});
list      = list(idx);
nfiles    = numel(list);
TAQmaster = cell(nfiles,1);
fmt       = ['%s %s %s ',repmat('%u8 ',1,9),'%s %f %f %s %u8 %f %*[^\n]'];
opts      = {'Delimiter',',','ReadVarNames',1,'CommentStyle',{'"','"'}};

% LOOP by file - 20 sec
for ii = 1:nfiles
    disp(list{ii})
    % Read in the whole master file
    TAQmaster{ii} = dataset('File',list{ii}, 'Format',fmt, opts{:});
    % Eventually rename DATEF to FDATE (in some files it changes)
    TAQmaster{ii}.Properties.VarNames = regexprep(TAQmaster{ii}.Properties.VarNames,'(?i)datef','FDATE');
    % Conversions
    TAQmaster{ii} = replacedata(TAQmaster{ii}, @logical, {'ETN','ETA','ETB','ETP','ETX','ETT','ETO','ETW','ITS'});
    TAQmaster{ii} = replacedata(TAQmaster{ii}, @int8, 'TYPE');
    TAQmaster{ii} = replacedata(TAQmaster{ii}, @uint32, {'FDATE','UOT'});
    % Convert to table
    TAQmaster{ii} = dataset2table(TAQmaster{ii});
end
delete(list{:})
%% Load .tab 
% Add .tab master files from dvd
fprintf('Unzipping .tab files\n')
list2     = unzip(fullfile(path2zips,tabtables),path2zips);
list2     = sort(list2);
fmt       = '%10c %30c %12c %c%c%c %*2c %c%c%c%c%c %c %4c %10c %4c %c %c %8c %*[^\n]';
names     = {'SYMBOL','NAME','CUSIP','ETN','ETA','ETB','ETP','ETX','ETT','ETO','ETW','ITS','ICODE','SHROUT','UOT','DENOM','TYPE','FDATE'};
nfiles2   = numel(list2);
TAQmaster = [TAQmaster; cell(nfiles2,1)];

% LOOP by file
for ii = 1:nfiles2
    disp(list2{ii})
    
    % Read in fixed width master files in a temporary variable
    fid = fopen(list2{ii});
    tmp = textscan(fid, fmt,'Whitespace','');
    fclose(fid);
    
    % Convert
    f = @(x,n) cellstr(x{n});
    g = @(x,n) num2cell(uint32(str2num(x{n})));
    h = @(x) [f(x,1) f(x,2) f(x,3)  num2cell(logical([x{4:12}]-'0'))  f(x,13) num2cell(str2num(x{14})) g(x,15) f(x,16) num2cell(uint8(x{17}-'0')) g(x,18)];
    TAQmaster{nfiles+ii} = cell2table(h(tmp),'VariableNames',names);
end
delete(list2{:})
%% Post-processing and saving

% Concatenate everything
TAQmaster = cat(1,TAQmaster{:});

% Unique records
TAQmaster = unique(TAQmaster,'stable');

% Save
save(fullfile(writeto, sprintf('%s_%s.mat',datestr(now,'yyyymmdd_HHMM'),'TAQmaster')), 'TAQmaster')
end

