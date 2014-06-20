function makeMasterFile(path2matfiles)
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
matlabpool close
toc

% Number of rows per block
len = cellfun(@(x) size(x,1),mst);

% Add number of file as 4th column to mst{2} and cat all mst
mst = arrayfun(@(x) [mst{x} dataset({repmat(x,len(x),1),'File'})],(1:numel(len))','un',0);
ids = cat(1,ids{:});
fprintf('Concatenating master records.\n')
tic
mst = cat(1,mst{:}); 
toc

% Re-index id to whole block
cumpos           = cumsum(len(1:end-1));
mst.Id           = [1; diff(mst.Id)~= 0];
mst.Id(cumpos+1) = 1;
mst.Id           = cumsum(mst.Id);

% Make unique tickers
[ids,~,long2unique] = unique(ids);
mst.Id              = long2unique(mst.Id);

% Sort according to id-date pair
mst = sortrows(mst,{'Id','Date'});

% Save
save(fullfile(path2matfiles,'master'),mstname, idsname,'-v6','-mat')
end