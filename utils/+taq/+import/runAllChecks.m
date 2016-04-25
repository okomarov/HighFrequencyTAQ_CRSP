function runAllChecks(matfolder,testfolder)
%% Day spillover
prob = checkDaySpillover(matfolder);
if ~isempty(prob)
    disp(prob)
    error('There are spillovers of daily data across multiple files.')
end

%% Master from-to continuity
dd         = dir(fullfile(matfolder,'*.idx'));
indexnames = {dd.name};
for ii = 1:numel(indexnames)
    load(fullfile(matfolder,indexnames{ii}),'-mat');
    if any(index.To(1:end-1)+1 ~= index.From(2:end))
        error('File %d has gaps between To and the next From.',ii)
    end
end

%% Compare price series
FMT    = '%s%u32%u8:%u8:%u8%f32%u32%u16%u16%s%c';
OPTS   = {'Delimiter',',','HeaderLines',1};
VNAMES = {'Symbol','Date','HH','MM','SS','Price','Size','G127','Corr','Cond','Ex'};

folders = dir(fullfile(testfolder));
folders = {folders.name};
folders = fullfile(testfolder, folders(~ismember(folders, {'.','..'})));
for f = folders
    zipfiles = dir(fullfile(f{1},'*.zip'));
    zipfiles = fullfile(f{1},{zipfiles.name});

    for z = zipfiles
        fname   = unzip(z{1},tempdir());
        fid     = fopen(fname{1});
        cleanup = onCleanup(@() myCleanup(fname{1},fid));

        manual = textscan(fid, FMT, OPTS{:});
        manual = table(manual{:},'VariableNames',VNAMES);

        automated = taq.getTrades('symbol',manual.Symbol{1},[manual.Date(1), manual.Date(end)],matfolder);
        time      = hhmmssmat2serial(automated.Time);
        idx       = in(time, hhmmssmat2serial([9 30 0; 16 0 0])');

        if ~isequal(manual.Price, automated.Price(idx))
            error('There are differences in the price series for "%s".',z{1})
        end
    end
end
end

function myCleanup(fname, fid)
fclose(fid);
delete(fname)
end
