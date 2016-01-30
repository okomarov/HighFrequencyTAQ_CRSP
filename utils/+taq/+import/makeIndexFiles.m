function makeIndexFiles(path2mst)
% MAKEINDEXFILES Groups all master records from single *.mat files into one master table

% Read .mat filenames
d = dir(fullfile(path2mst,'*.mst'));

% Preallocate
nfiles = numel(d);

[mstSymb, mstDate] = preallocate(nfiles);

for f = 1:nfiles
    disp(f/nfiles*100)
    s = load(fullfile(path2mst,d(f).name),'-mat');

    % Symbol-file index
    mstSymb = mapList(mstSymb, s.ids, f);

    % Date-file index
    mstDate = mapList(mstDate, unique(s.mst.Date), f);

end

% Drop excess pre-allocation
mstSymb = dropEmptyPreallocated(mstSymb);
mstDate = dropEmptyPreallocated(mstDate);

% Organize in sorted tables
mstSymb = table(mstSymb.Universe, mstSymb.FileList, 'VariableNames',{'Symbol','File'});
mstDate = table(mstDate.Universe, mstDate.FileList, 'VariableNames',{'Date','File'});
mstSymb = sortrows(mstSymb,'Symbol');
mstDate = sortrows(mstDate,'Date');

save(fullfile(path2mst,'master'),'mstSymb', 'mstDate','-v6','-mat')
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
    numNew           = nnz(~imember);
    n                = getUniverseSize(s.Universe);
    newPosInUni      = n+1:n+numNew;
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
fileClass = 'uint16';
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
