function makeMasterFile(path2matfiles)
% MAKEMASTERFILE Groups all master records from single *.mat files into one master table

mstname = 'mst';
idsname = 'ids';
% Read .mat filenames
d = dir(fullfile(path2matfiles,'*.mat'));

% Open matlabpool
if matlabpool('size') == 0
    matlabpool('open', 4, 'AttachedFiles',{'.\utils\poolStartup.m'})
end

% Preallocate
nfiles    = numel(d);
[mst,ids] = deal(cell(nfiles,1));
tic
% Load mst and ids
fprintf('Loading .mat files.\n')
parfor f = 1:nfiles
    disp(f)
    s = load(fullfile(path2matfiles,d(f).name),mstname, idsname);
    mst{f}  = s.(mstname);
    ids{f}  = s.(idsname);
end
toc

% Number of rows per block
lenmst = cellfun(@(x) size(x,1),mst);
lenids = cumsum(cellfun(@(x) size(x,1),ids));

% Add number of file as 4th column to mst{2} and cat all mst
mst = arrayfun(@(x) [mst{x} dataset({repmat(x,lenmst(x),1),'File'})],(1:nfiles)','un',0);

% Re-index Id
parfor f = 2:nfiles
    mst{f}.Id = mst{f}.Id + lenids(f-1);
end
matlabpool close

% Concatenate
ids = cat(1,ids{:});
fprintf('Concatenating master records.\n')
tic
mst = cat(1,mst{:}); 
toc

% Make unique tickers
[ids,~,long2unique] = unique(ids);
mst.Id              = long2unique(mst.Id);

% Drop tickers without data
ikeep = ismember(1:max(mst.Id), mst.Id);
ids   = ids(ikeep);

% Sort according to id-date pair
mst = sortrows(mst,{'Id','Date'});

% Save
save(fullfile(path2matfiles,'master'),mstname, idsname,'-v6','-mat')
end