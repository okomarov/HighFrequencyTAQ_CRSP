function matnum = splitMatTrades(path2data, outdir, nrecords, matnum, opt)
% SPLITMATTRADES Splits mat files by number of records
%
%   SPLITMATTRADES(PATH2DATA, [OUTDIR], [NRECORDS], [MATNUM], [OPT])
%       - PATH2DATA directory where .mat and .mst files to split are
%       - OUTDIR defaults to PATH2DATA if unspecified
%       - NRECORDS defaults to 5e6
%       - MATNUM saved mat/mst files will be numbered from MATNUM+1 (default: 0)

if nargin < 2 || isempty(outdir)
    outdir = path2data;
end
if nargin < 3 || isempty(nrecords)
    nrecords = 5e6;
end
if nargin < 4, matnum = 0; end
if nargin < 5
    opt.Fmt  = 'T%05d';
end

% Read .mat filenames
list_mat = dir(fullfile(path2data,'*.mat'));
list_mst = dir(fullfile(path2data,'*.mst'));

% Preallocate
[resdata, resmst,resids] = deal([]);
nfiles                   = numel(list_mat);
for f = 1:nfiles
    disp(list_mat(f).name)
    % Load
    sd = load(fullfile(path2data,list_mat(f).name));
    sm = load(fullfile(path2data,list_mst(f).name),'-mat');

    [sd.data,sm.mst,sm.ids] = prependResidualData(sd.data,sm.mst,sm.ids, resdata,resmst,resids);

    % Bin chunks of data
    nrows = size(sd.data,1);
    edges = 0:nrecords:nrows;
    if edges(end) < nrows
        edges = [edges, nrows];
    end
    [~,~,bin] = histcounts(sm.mst.To,edges);

    % Save each chunk into its own file
    for ii = 1:numel(edges) - 2
        matnum = matnum + 1;

        % In case of crash
        save(fullfile(outdir,'bck.mat'),'edges','ii','f','matnum')

        idx            = bin == ii;
        [data,mst,ids] = getDataMstIds(sd.data, sm.mst(idx,:), sm.ids);

        datafname = fullfile(outdir, sprintf([opt.Fmt, '.mat'],matnum));
        mstfname  = fullfile(outdir, sprintf([opt.Fmt, '.mst'],matnum));
        save(datafname,'data','-v7.3')
        save(mstfname ,'mst','ids','-v6')
    end
    % Carry over residual records
    idx                     = bin == ii+1;
    [resdata,resmst,resids] = getDataMstIds(sd.data, sm.mst(idx,:), sm.ids);
    save(fullfile(outdir,'rsd.mat'),'resdata','resmst','resids')
end
% Last residual records get their own
if isempty(resdata)
    return
end
matnum    = matnum + 1;
data      = resdata;
mst       = resmst;
ids       = resids;
datafname = fullfile(outdir, sprintf([opt.Fmt, '.mat'],matnum));
mstfname  = fullfile(outdir, sprintf([opt.Fmt, '.mst'],matnum));
save(datafname,'data','-v7.3')
save(mstfname ,'mst','ids','-v6')
end

function [data,mst,ids] = prependResidualData(data,mst,ids, resdata, resmst,resids)
if isempty(resdata)
    return
end
% Shift from/to
nrows    = size(resdata,1);
mst.From = mst.From + nrows;
mst.To   = mst.To + nrows;

% Remap ids
ids_new      = unique([resids; ids]);
[~,id_resid] = ismember(resids(resmst.Id), ids_new);
resmst.Id    = id_resid;
[~,id_ids]   = ismember(ids(mst.Id), ids_new);
mst.Id       = id_ids;

ids  = ids_new;
data = [resdata; data];
mst  = [resmst; mst];
end

function [data,mst,ids] = getDataMstIds(data, mst, ids)
if isempty(mst)
    [data, mst, ids] = deal([]);
    return
end
if ~issorted(mst.From)
    error('mst should be sorted');
end
[pos,~,id] = unique(mst.Id);
ids        = ids(pos);
mst.Id     = id;
data       = data(mst.From(1):mst.To(end),:);
mst.To     = mst.To - mst.From(1) + 1;
mst.From   = mst.From - mst.From(1) + 1;
end