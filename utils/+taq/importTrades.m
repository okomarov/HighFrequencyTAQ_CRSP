function importTrades(rootfolder)
%   IMPORTTRADES Imports TAQ raw trades located in the raw folder
%
% Assumes the following folder structure:
%
%     rootfolder/
%     +-- manual/
%     |     199301_01.zip
%     |     ...
%     |     201005_03.zip
%     |
%     +-- dvd/
%     |   +-- 201006/
%     |   |     CDA.zip
%     |   |     CDB.zip
%     |   |     ...
%     |   +-- .../
%     |   \-- 201012/
%     |         CDA.zip
%     |         CDB.zip
%     |         ...
%     |
%     \-- automated/
%           CT_20110103.zip
%             ...
%
% OUTPUT: 
% The .mat files will contain the data records (table), while the .idx
% files will contain the index records (table) and the unique list of tickers
% (cellstring) that maps directly to the numeric ID.
%
% - data:
%       .Time [hh, mm, ss]                          uint8
%       .Price                                      single
%       .Volume                                     uint32
%       .G127_Correction [G127 rule,  correction]   uint16
%       .Condition                                  2char
%       .Exchange                                   1char
%
% - symbol, unique list of tickers                     cellstring
%
% - index, master of id-date pairs mapping to data and ids:
%       .Id, map to ids                             uint16
%       .Date, yyyymmdd                             uint32
%       .From, starting row in data                 uint32
%       .To,   ending row in data                   uint32
%
% For details on the data and its fields see <a href="http://www.nyxdata.com/Data-Products/Monthly-TAQ-DVD">Monthly TAQ DVD</a>.
recycle off
diary(fullfile(rootfolder,'log.txt'))

% Create rootfolder/mat if does not exist
matfolder = fullfile(rootfolder,'..','mat');
if ~isdir(matfolder)
    mkdir(matfolder);
end

%% Imports

% Import manual csv
% =================
% Note: 7388 files at ~5e6 records per file
folder = fullfile(rootfolder, 'manual');
matnum = taq.import.tradesCSV(folder,matfolder,true);

% Import dvd 
% ==========
% NOTE: 146 files, and 791 after splitting into 5e6 records per file

% Data spillovers across files does not happen so import as they come and
% chunk into smaller pieces later
folder = fullfile(rootfolder, 'dvd');
% Import into tmp folder
tmpmat = fullfile(matfolder, 'tmp');
if ~isdir(tmpmat)
    mkdir(tmpmat);
end
taq.import.tradesDVD(folder,tmpmat,matnum);
% Chunk into smaller .mat files
matnum = taq.import.splitMatTrades(tmpmat,matfolder,[],matnum);

% Import automated csv
% ====================
folder = fullfile(rootfolder, 'automated');
taq.import.tradesCSV(folder,matfolder,false, matnum);

%% Checks
testfolder = fullfile(rootfolder,'test');
taq.import.runAllChecks(matfolder,testfolder);

% Tests against custom excels (need TAQ.getTrades)
diary off
end