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
folder = fullfile(rootfolder, 'manual');
importTradesCSV(folder,matfolder,1);

% Import dvd 
folder = fullfile(rootfolder, 'dvd');
importTradesDVD(folder,matfolder);

% Import automated csv
folder = fullfile(rootfolder, 'automated');
importTradesCSV(folder,matfolder,false);

diary off
end