function makeMasterFile(path2matfiles, isdebug)
% MAKEMASTERFILE Groups all master records from single *.mat files into one master table

if nargin < 2 || isempty(isdebug), isdebug = false; end

mstname = 'mst';
idsname = 'ids';
% Read .mat filenames
d = dir(fullfile(path2matfiles,'*.mat'));

fpath = fullfile(fileparts(mfilename('fullpath')), 'poolStartup.m');
poolStartup(4, 'AttachedFiles',{fpath},'debug',isdebug)

% Preallocate
nfiles    = numel(d);
[mst,ids] = deal(cell(nfiles,1));
tic
% Load mst and ids
fprintf('Loading .mat files.\n')
parfor f = 1:nfiles
    disp(f)
    s      = load(fullfile(path2matfiles,d(f).name),mstname, idsname);
    if isa(s.(mstname), 'dataset')
        mst{f} = dataset2table(s.(mstname));
    else
        mst{f} = s.(mstname);
    end
    ids{f} = s.(idsname);
end
delete(gcp)
toc

% Match id to unique ids
% Note: keep the matching algo simple (had enough of debugging)
lenids      = cellfun('prodofsize', ids);
[ids,~,pos] = unique(cat(1,ids{:}));
pos         = uint16(pos);
pos         = mat2cell(pos,lenids);
for f = uint16(1:nfiles)
   mst{f}.Id   = pos{f}(mst{f}.Id);
   mst{f}.File = repmat(f, size(mst{f},1),1);
end

% Concatenate all
fprintf('Concatenating records.\n')
mst = cat(1,mst{:}); 

% Drop tickers without data
idrop = ~ismember(1:numel(ids), mst.Id);
ids(idrop) = {''};

% Convert From/To to uint32
mst.From = uint32(mst.From);
mst.To   = uint32(mst.To);

% Sort according to id-date pair
mst = sortrows(mst,{'Id','Date'});

% Save
save(fullfile(path2matfiles,'master'),mstname, idsname,'-v6','-mat')
end