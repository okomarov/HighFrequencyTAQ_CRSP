function matnum = splitMatTrades(path2data, outdir, nrecords, matnum, opt)
% SPLITMATTRADES Splits mat files by number of records
%
%   SPLITMATTRADES(PATH2DATA, [OUTDIR], [NRECORDS], [MATNUM], [OPT])
%       - PATH2DATA directory where .mat and .idx files to split are
%       - OUTDIR defaults to PATH2DATA if unspecified
%       - NRECORDS defaults to 5e6
%       - MATNUM saved mat/idx files will be numbered from MATNUM+1 (default: 0)

if nargin < 2 || isempty(outdir)
    outdir = path2data;
end
if nargin < 3 || isempty(nrecords)
    nrecords = 5e6;
end
if nargin < 4, matnum = 0; end
if nargin < 5
    opt.Fmt = 'T%05d';
end

% Read .mat filenames
list_mat = dir(fullfile(path2data,'*.mat'));
list_idx = dir(fullfile(path2data,'*.idx'));

% Preallocate
[resdata, resindex,ressymbol] = deal([]);
nfiles                        = numel(list_mat);
for f = 1:nfiles
    disp(list_mat(f).name)
    % Load
    sd = load(fullfile(path2data,list_mat(f).name));
    sm = load(fullfile(path2data,list_idx(f).name),'-mat');

    [sd.data,sm.index,sm.symbol] = prependResidualData(sd.data,sm.index,sm.symbol, resdata,resindex,ressymbol);

    % Bin chunks of data
    nrows = size(sd.data,1);
    edges = 0:nrecords:nrows;
    if edges(end) < nrows
        edges = [edges, nrows];
    end
    [~,~,bin] = histcounts(sm.index.To,edges);

    % Save each chunk into its own file
    for ii = 1:numel(edges) - 2
        matnum = matnum + 1;

        % In case of crash
        save(fullfile(outdir,'bck.mat'),'edges','ii','f','matnum')

        idx                 = bin == ii;
        [data,index,symbol] = getDataMstIds(sd.data, sm.index(idx,:), sm.symbol);

        datafname  = fullfile(outdir, sprintf([opt.Fmt, '.mat'],matnum));
        indexfname = fullfile(outdir, sprintf([opt.Fmt, '.idx'],matnum));
        save(datafname,'data','-v7.3')
        save(indexfname ,'index','symbol','-v6')
    end
    % Carry over residual records
    idx                          = bin == ii+1;
    [resdata,resindex,ressymbol] = getDataMstIds(sd.data, sm.index(idx,:), sm.symbol);
    save(fullfile(outdir,'rsd.mat'),'resdata','resindex','ressymbol')
end
% Last residual records get their own
if isempty(resdata)
    return
end
matnum     = matnum + 1;
data       = resdata;
index      = resindex;
symbol     = ressymbol;
datafname  = fullfile(outdir, sprintf([opt.Fmt, '.mat'],matnum));
indexfname = fullfile(outdir, sprintf([opt.Fmt, '.idx'],matnum));
save(datafname,'data','-v7.3')
save(indexfname ,'index','symbol','-v6')
end

function [data,index,symbol] = prependResidualData(data,index,symbol, resdata, resindex,ressymbol)
if isempty(resdata)
    return
end
% Shift from/to
nrows      = size(resdata,1);
index.From = index.From + nrows;
index.To   = index.To + nrows;

% Remap symbols
symbol_new    = unique([ressymbol; symbol]);
[~,id_resid]  = ismember(ressymbol(resindex.Id), symbol_new);
resindex.Id   = id_resid;
[~,id_symbol] = ismember(symbol(index.Id), symbol_new);
index.Id      = id_symbol;

symbol = symbol_new;
data   = [resdata; data];
index  = [resindex; index];
end

function [data,index,symbol] = getDataMstIds(data, index, symbol)
if isempty(index)
    [data, index, symbol] = deal([]);
    return
end
if ~issorted(index.From)
    error('Index table should be sorted');
end
[pos,~,id] = unique(index.Id);
symbol     = symbol(pos);
index.Id   = id;
data       = data(index.From(1):index.To(end),:);
index.To   = index.To - index.From(1) + 1;
index.From = index.From - index.From(1) + 1;
end