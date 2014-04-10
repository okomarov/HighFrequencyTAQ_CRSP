function data = loadPriceTAQ(ticker,from,to)

if nargin < 2 || isempty(from), from = 0; end
if nargin < 3 || isempty(to),   to = inf; end

master = load(fullfile(cd, '.\data\TAQ\master'),'-mat');
id = find(strcmpi(master.ids,ticker));
if isempty(id)
    warning('ticker ''%s'' unmatched.',upper(ticker))
    return
end

master = master.mst(master.mst.Id == id,:);
master = master(master.Date <= to & master.Date >= from,:);
blocks = master.To - master.From + 1;

data   = zeros(sum(blocks),10);

h = waitbar(0, sprintf('Loading ticker ''%s''', upper(ticker)));
files  = unique(master.File);
nfiles = numel(files);
c = 0;

s = load(fullfile(cd, '.\data\TAQ', sprintf('T%04d.mat',files(1))),'data');
for ii = 1:size(master,1)-1
        
    idx = master.File == files(ii);
    for jj = 1:numel(master.File(idx))
        c = c + blocks()
    end
    if master.File(ii) < master.File(ii+1)
        s = load(fullfile(cd, '.\data\TAQ', sprintf('T%04d.mat',files(ii+1))),'data');
        waitbar(ii/nfiles)
end
close(h)
end