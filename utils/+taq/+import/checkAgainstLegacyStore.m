function checkAgainstLegacyStore(path2legacy, matfolder)
matfolder   = 'I:\Box Sync\TAQ\mat';
path2legacy = 'data\TAQ\master';
master      = load(path2legacy, '-mat');

% Check old master versus new
dd     = dir(fullfile(matfolder,'*.idx'));
nfiles = numel(dd);
new    = cell(nfiles,1);
for ii = 1:nfiles
    tmp          = load(fullfile(matfolder,dd(ii).name),'-mat');
    [~,pos]      = ismember(tmp.symbol, master.ids);
    tmp.index.Id = pos(tmp.index.Id);
    new{ii}      = tmp.index(:,{'Id','Date'});
end
new = cat(1,new{:});
new = sortrows(new,{'Id','Date'});
if ~isequal(new, master.mst(:,{'Id','Date'}))
    error('The legacy master is different from the newly created master.')
end

%% Check random series
NUM_TESTS = 10000;

for ii = 1:NUM_TESTS
    obs    = datasample(master.mst,1);
    symbol = master.ids{obs.Id};
    new    = taq.getTrades('symbol',symbol,[obs.Date, obs.Date],matfolder,false);
    old    = getTaqData('symbol',symbol,obs.Date, obs.Date);
    if isempty(new) || ~isequal(new.Price, old.Price)
        error('There legacy price series is different from the newly retrieved price series.');
    end
end
end
