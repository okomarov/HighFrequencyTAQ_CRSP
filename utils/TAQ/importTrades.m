function importTrades(rootfolder)
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
% The .mat files will contain the data records (table), while the .mst
% files will contain the mst records (table) and the unique list of tickers
% (cellstring) that maps to mst.Id.
%
% - data:
%       .Time [hh, mm, ss]                          uint8
%       .Price                                      single
%       .Volume                                     uint32
%       .G127_Correction [G127 rule,  correction]   uint16
%       .Condition                                  2char
%       .Exchange                                   1char
%
% - ids, unique list of tickers                     cellstring
%
% - mst, master of id-date pairs mapping to data and ids:
%       .Id, map to ids                             uint16
%       .Date, yyyymmdd                             uint32
%       .From, starting row in data                 uint32
%       .To,   ending row in data                   uint32
%
% For details on the data and its fields see <a href="http://www.nyxdata.com/Data-Products/Monthly-TAQ-DVD">Monthly TAQ DVD</a>.
recycle off
diary(fullfile(rootfolder,'log.txt'))

% Create rootfolder/mat if does not exist
matfolder = fullfile(rootfolder,'mat');
if ~isdir(matfolder)
    mkdir(matfolder);
end

% Import manual csv
% =================
% Note: 7388 files at ~5e6 records per file
folder = fullfile(rootfolder, 'manual');
matnum = import.tradesCSV(folder,matfolder,true);

% Import dvd 
% ==========
% NOTE: 146 files, and 791 after splitting into 5e6 records per file

% Data spillovers across files does not happen so import as they come and
% chunk into smaller pieces later
folder = fullfile(rootfolder, 'dvd');
% Import into tmp folder
tmpmat = fullfile(rootfolder, 'tmp');
if ~isdir(tmpmat)
    mkdir(tmpmat);
end
import.tradesDVD(folder,tmpmat,matnum);
% Chunk into smaller .mat files
matnum = import.splitMatTrades(tmpmat,fullfile(tmpmat,'split'),[],matnum);

% Import automated csv
% ====================
folder = fullfile(rootfolder, 'automated');
import.tradesCSV(folder,matfolder,false, matnum);

diary off
end