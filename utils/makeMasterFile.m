function makeMasterFile(path2matfiles, isdebug)
% MAKEMASTERFILE Groups all master records from single *.mat files into one master table

if nargin < 2 || isempty(isdebug), isdebug = false; end

mstname = 'mst';
idsname = 'ids';
% Read .mat filenames
d = dir(fullfile(path2matfiles,'*.mat'));

% Open matlabpool
if isempty(gcp('nocreate')) && ~isdebug
    parpool(4, 'AttachedFiles',{'.\utils\poolStartup.m'})
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
    mst{f}     = sortrows(s.(mstname),'Id');
    mst{f}.Id  = int64(mst{f}.Id);
    ids{f}     = s.(idsname);
end
delete(gcp)
toc

% Number of rows per block
lenmst = cellfun(@(x) size(x,1),mst);

% Shift in Ids from a block to another, i.e. not all ids are in mst.Id
lenids    = cellfun(@(x) size(x,1),ids);
maxids    = cellfun(@(x) max(x.Id),mst);
endshifts = cast(lenids,'like',maxids) - maxids + 1;

% Concatenate all
fprintf('Concatenating records.\n')
ids = cat(1,ids{:});
mst = cat(1,mst{:}); 

% Add number of file to mst
mst.File(:,1) = RunLength(uint16(1:nfiles), lenmst);

% Re-index id to whole block
cumpos           = cumsum(lenmst(1:end-1));
mst.Id           = [1; diff(mst.Id)];
mst.Id(cumpos+1) = endshifts(1:end-1);
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