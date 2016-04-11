function matnum = tradesDVD(path2main,outdir, matnum, opt)
% TRADESDVD Import zipped .csv of TAQ trades provided by NYSE through DVDs
%
%   TRADESCSV(PATH2MAIN)
%       PATH2MAIN points to the main folder with the subfolders containing
%       zipped .CSVs (see importTrades for the folder structure)
%
% For detailed help, see tradesCSV (when applicable), e.g. OPT has only the
% Fmt field.
%
% See also: TRADESCSV, IMPORTTRADES, TEXTSCAN
if nargin < 3, matnum = 0; end

if nargin < 4
    opt.Fmt = 'T%05d';
end

% Get monthly subfolders with zipped data in CDA.zip, CDB.zip, ...
d          = dir(path2main);
d          = d([d.isdir]);
subfolders = {d.name};
subfolders = subfolders(~ismember(subfolders,{'.','..'}));

% LOOP by subfolder
for ii = 1:numel(subfolders)
    datadir  = fullfile(path2main,subfolders{ii});
    s        = dir(fullfile(datadir,'*.zip'));
    zipnames = sort({s.name});

    % LOOP by zip file
    for jj = 1:numel(zipnames)
        matnum  = matnum+1;
        zipfile = fullfile(datadir,zipnames{jj});
        disp(zipfile)
        try
            % Unzip trades, e.g. T201006A
            zipcontents   = char(getZipEntries(zipfile));
            iTradeFiles   = zipcontents(:,1) == 'T';
            files2extract = cellstr(zipcontents(iTradeFiles,:));
            extracted     = unzipfiles(zipfile, files2extract, tempdir());
            extracted     = fullfile(tempdir(), extracted);
            cleanup       = onCleanup(@() delete(extracted{:}));

            % Import index records, extension .IDX
            idx            = ~cellfun('isempty',regexpi(extracted, '\.idx$'));
            [symbol,index] = importIDX(extracted{idx});

            % Read in data, extension .BIN
            data = importBIN(extracted{~idx});

            % 3. Save data and cleanup
            datafname  = fullfile(outdir, sprintf([opt.Fmt '.mat'],matnum));
            indexfname = fullfile(outdir, sprintf([opt.Fmt '.idx'],matnum));
            save(datafname,'data','-v7.3')
            save(indexfname ,'index','symbol','-v6')
            fprintf('%-40s%s\n',sprintf([opt.Fmt '.mat'],matnum),datestr(now,'dd HH:MM:SS'))
            delete(cleanup)

        catch err
            filename = fullfile(path2main, sprintf('%s_err.mat',datestr(now,'yyyymmdd_HHMM')));
            save(filename,'err')
            rethrow(err)
        end
    end
end
% Chunk mat files
end

function entries = getZipEntries(zipfilename)
jZipFile = java.util.zip.ZipFile(zipfilename);
jEntries = jZipFile.entries;
entries  = {};
c        = 0;
while jEntries.hasMoreElements
    c          = c + 1;
    entries{c} = char(jEntries.nextElement);
end
end

function [symbol,index] = importIDX(filename)
% symbol | date | from | to
%   10      4      4      4
fid     = fopen(filename);
cleanup = onCleanup(@() fclose(fid));

% Retrieve the number of rows
fseek(fid,0,'eof');
cols  = 22;
nrows = ftell(fid)/cols;

% Read in whole data
fseek(fid,0,'bof');
tmp = fread(fid,[cols,nrows],'*uint8')';

% Get symbols and mapping (no sorting)
[symbol, ~, id] = unique(cellstr(char(tmp(:,1:10))),'stable');

% Create index table
fun32 = @(x,pos) typecast(reshape(x(:,pos)',[],1),'uint32');
index = table(uint16(id)      ,...
              fun32(tmp,11:14),...
              fun32(tmp,15:18),...
              fun32(tmp,19:22),  'VariableNames',{'Id','Date','From','To'});
end

function data = importBIN(filename)
% Seconds | Price | Volume | G127 and Correction | Condition | Exchange
%    4        4       4              4                 2          1

fid     = fopen(filename);
cleanup = onCleanup(@() fclose(fid));

% Retrieve the number of rows
fseek(fid,0,'eof');
cols  = 19;
nrows = ftell(fid)/cols;

% Read in whole data
fseek(fid,0,'bof');
tmp = fread(fid,[19,nrows],'*uint8')';

% Create data table
timefun = @(x) uint8(fix([mod(x,86400)/3600, mod(x,3600)/60, mod(x,60)]));
fun16   = @(x,pos) typecast(reshape(x(:,pos)',[],1),'uint16');
fun32   = @(x,pos) typecast(reshape(x(:,pos)',[],1),'uint32');
data    = table(timefun(double(fun32(tmp,1:4))),...
                single(fun32(tmp,5:8))/1e4     ,...
                fun32(tmp,9:12)                ,...
                reshape(fun16(tmp,13:16),2,[])',...
                char(tmp(:,17:18))             ,...
                char(tmp(:,   19)), 'VariableNames', {'Time', 'Price','Volume','G127_Correction','Condition','Exchange'});
end