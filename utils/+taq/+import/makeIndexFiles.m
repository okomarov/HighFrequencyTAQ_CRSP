function makeIndexFiles(path2idx)
% MAKEINDEXFILES Groups all index records from single *.idx files into one index table

PREFIX_INDEX = 'index_';

% Read .idx filenames
d = dir(fullfile(path2idx,'*.idx'));

% Preallocate
nfiles = numel(d);

[indexSymb, indexDate] = preallocate(nfiles);

for f = 1:nfiles
    disp(f/nfiles*100)
    s = load(fullfile(path2idx,d(f).name),'-mat');

    % Symbol-file index
    indexSymb = mapList(indexSymb, s.symbol, f);

    % Date-file index
    indexDate = mapList(indexDate, unique(s.index.Date), f);
end

% Drop excess pre-allocation
indexSymb = dropEmptyPreallocated(indexSymb);
indexDate = dropEmptyPreallocated(indexDate);

% Organize in containers
symbol = containers.Map(upper(indexSymb.Universe), indexSymb.FileList);
id     = containers.Map(uint32(1:numel(indexSymb.Universe)), indexSymb.FileList);
date   = containers.Map(indexDate.Universe, indexDate.FileList);

save(fullfile(path2idx,[PREFIX_INDEX, 'symbol']),'symbol','-v6','-mat')
save(fullfile(path2idx,[PREFIX_INDEX, 'id'])    ,'id'    ,'-v6','-mat')
save(fullfile(path2idx,[PREFIX_INDEX, 'date'])  ,'date'  ,'-v6','-mat')
end

function s = mapList(s, records, f)
% Adds the file number to the file list of a member of the universe.
%
%   MAPLIST(S, RECORDS, F)
%       S is a structure with:
%           .Universe - cell array of symbols or uints (e.g. dates)
%           .FileList - array of lists with the file numbers
%           .FirstEmptyInList - array of positions pointing, for each list,
%               to the first empty spot where we can add the file number
%       RECORDS has the members to which I will add the file number F
%       F is the file number containing RECORDS
%
%   If a RECORD is not yet a member, grows the UNIVERSE. For performance,
%   the universe, the list and array of positions are all pre-allocated.

[imember,existingPos] = ismember(records,s.Universe);

% Existing members
if any(imember)
    existingPos = existingPos(imember);
    for ii = 1:nnz(imember)
        p = existingPos(ii);
        s.FileList{p}(s.FirstEmptyInList(p)) = f;
    end
    % Grow positions
    s.FirstEmptyInList(existingPos) = s.FirstEmptyInList(existingPos)+1;
end

% Add new members to universe and their lists
if any(~imember)
    numNew                  = nnz(~imember);
    n                       = getUniverseSize(s.Universe);
    newPosInUni             = n+1:n+numNew;
    s.Universe(newPosInUni) = records(~imember);
    for ii = 1:numNew
        p = newPosInUni(ii);
        s.FileList{p}(s.FirstEmptyInList(p)) = f;
    end
    % Add positions to new members
    s.FirstEmptyInList(newPosInUni) = s.FirstEmptyInList(newPosInUni)+1;
end
end

function [s1,s2] = preallocate(nfiles)
fileClass           = 'uint16';
% Symbol
N                   = 4e4;
s1.Universe         = repmat({''},N,1);
s1.FileList         = repmat({zeros(1,nfiles,fileClass)},N,1);
s1.FirstEmptyInList = ones(N,1,fileClass);
% Date
N                   = 1e4;
s2.Universe         = zeros(N,1,'uint32');
s2.FileList         = repmat({zeros(1,nfiles,fileClass)},N,1);
s2.FirstEmptyInList = ones(N,1,fileClass);
end

function s = dropEmptyPreallocated(s)
if iscell(s.Universe)
    ikeep = ~cellfun('isempty',s.Universe);
else
    ikeep = s.Universe ~= 0;
end
s.Universe = s.Universe(ikeep);
s.FileList = cellfun(@(x) x(x~=0), s.FileList(ikeep),'un',0);
end

function n = getUniverseSize(universe)
% Gets size of a preallocated UNIVERSE array

% Position of first empty (preallocated) element - 1
if iscell(universe)
    n = find(cellfun('isempty',universe),1,'first')-1;
else
    n = find(universe == 0,1,'first')-1;
end
end
