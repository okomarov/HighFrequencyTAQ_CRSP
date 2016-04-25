function addPermnoToIndexFiles(path2idx, path2mst)
if nargin < 2
    rootfld  = fileparts(path2idx);
    path2mst = fullfile(rootfld,'master');
end

% Read .idx filenames
d = dir(fullfile(path2idx,'*.idx'));

nfiles = numel(d);
for f = 300:nfiles
    disp(f/nfiles*100)

    fname = fullfile(path2idx,d(f).name);
    s     = load(fname,'-mat');

    cusip          = getAllCusip(s, path2mst);
    s.index.Permno = crsp.getPermno(cusip);

    save(fname,'-struct','s','-mat','-v6')
end
end

function cusip = getAllCusip(s, path2mst)
% Get symbol according to numeric id
first_id = [true; logical(diff(uint64(s.index.Id)))];
symbols  = s.symbol(s.index.Id(first_id));

% Cache dates
dates = cache2cell(s.index.Date, cumsum(first_id));

N     = numel(symbols);
cusip = cell(N,1);
parfor ii = 1:N
    try
        cusip{ii} = taq.getCusip(symbols{ii}, dates{ii}, path2mst);
    catch ME
        if ~strcmpi(ME.identifier,'MATLAB:load:couldNotReadFile')
            disp(ME.message)
            rethrow(ME)
        end
    end
    if isempty(cusip{ii})
        cusip{ii} = repmat(' ',numel(dates{ii}),8);
    end
end
cusip = cat(1,cusip{:});
end