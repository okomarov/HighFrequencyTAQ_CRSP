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
    mst{f}  = sortrows(s.(mstname),'Id');
    mst{f}  = int64(mst{f}.Id);
    ids{f}  = s.(idsname);
end
matlabpool close
toc

% Number of rows per block
lenmst = cellfun(@(x) size(x,1),mst);

% Add number of file as 4th column to mst{2} and cat all mst
mst = arrayfun(@(x) [mst{x} dataset({repmat(x,lenmst(x),1),'File'})],uint16((1:nfiles))','un',0);
ids = cat(1,ids{:});
fprintf('Concatenating master records.\n')
tic
mst = cat(1,mst{:}); 
toc

% Re-index id to whole block
cumpos           = cumsum(lenmst(1:end-1));
mst.Id           = [1; diff(mst.Id)];
mst.Id(cumpos+1) = 1;
mst.Id           = cumsum(mst.Id);

% Make unique tickers
[ids,~,long2unique] = unique(ids);
mst.Id              = uint16(long2unique(mst.Id));

% Drop tickers without data
idrop = ~ismember(1:numel(ids), mst.Id);
ids(idrop) = {''};

% Sort according to id-date pair
mst = sortrows(mst,{'Id','Date'});

% Save
save(fullfile(path2matfiles,'master'),mstname, idsname,'-v6','-mat')
end