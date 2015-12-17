function makeIndexFiles(path2mst)
% MAKEINDEXFILES Groups all master records from single *.mat files into one master table

% Read .mat filenames
d = dir(fullfile(path2mst,'*.mst'));

% Preallocate
nfiles = numel(d);

N            = 4e4;
Symbol       = repmat({''},N,1);
fileList1    = repmat({NaN(1,nfiles)},N,1);
posLastFile1 = ones(N,1);

N            = 1e4;
Date         = zeros(N,1,'uint32');
fileList2    = repmat({NaN(1,nfiles)},N,1);
posLastFile2 = ones(N,1);

for f = 1:nfiles
    disp(f)
    s = load(fullfile(path2mst,d(f).name),'-mat');

    % Symbol-file index
    [Symbol, fileList1, posLastFile1] = mapList(Symbol, s.ids, fileList1, posLastFile1, f);

    % Date-file index
    [Date, fileList2, posLastFile2] = mapList(Date, unique(s.mst.Date), fileList2, posLastFile2, f);
    
end

% Drop excess pre-allocation
ikeep     = ~cellfun('isempty',Symbol);
Symbol    = Symbol(ikeep);
fileList1 = cellfun(@(x) x(~isnan(x)), fileList1(ikeep),'un',0);

ikeep     = Date ~= 0;
Date      = Date(ikeep);
fileList2 = cellfun(@(x) x(~isnan(x)), fileList2(ikeep),'un',0);

% Organize in sorted tables
Symbols = table(Symbol(ikeep), fileList1(ikeep), 'VariableNames',{'Symbol','File'});
Dates   = table(Date(ikeep), fileList2(ikeep), 'VariableNames',{'Date','File'});
Symbols = sortrows(Symbols,'Symbol');
Dates   = sortrows(Dates,'Date');

save(fullfile(path2mst,'master'),mstname, idsname,'-v6','-mat')
end

function [universe, attribute, lasPosAttr] = mapList(universe, newset, attribute, lasPosAttr, f)
[idx,pos] = ismember(newset,universe);

% Grow list of existing attributes
if any(idx)
    pos = pos(idx);
    for ii = 1:numel(pos)
        p                           = pos(ii);
        attribute{p}(lasPosAttr(p)) = f;
    end
    % Grow positions
    lasPosAttr(pos) = lasPosAttr(pos)+1;
end

% Grow universe and add attributes
if any(~idx)
    numNew = nnz(~idx);
    if iscell(universe)
        uniSize = find(cellfun('isempty',universe),1,'first')-1;
    else
        uniSize = find(universe==0,1,'first')-1;
    end
       
    pos           = uniSize+1:uniSize+numNew;
    universe(pos) = newset(~idx);
    for ii = 1:numNew
        p                           = pos(ii);
        attribute{p}(lasPosAttr(p)) = f;
    end
    lasPosAttr(pos) = lasPosAttr(pos)+1;
end
end