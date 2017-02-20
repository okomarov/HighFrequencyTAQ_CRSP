function res = getBadPricesCount(path2data, outname, cleanOpts, iterOpts)
% Parse options
if nargin < 3 || isempty(cleanOpts)
    cleanOpts = struct('ExcludeBadDays',false,'BadPriceMultiplier', 10,'TerminateEarly',true,'ConsolidateTimestampType',[]);
else
    p              = inputParser();
    p.StructExpand = true;
    addParameter(p,'BadPriceMultiplier',10);
    addParameter(p,'ExcludeBadDays',false);
    addParameter(p,'TerminateEarly',true);
    addParameter(p,'ConsolidateTimestampType',[]);
    p.parse(cleanOpts);
    cleanOpts      = p.Results;
end

% Parse iterOpts
if nargin < 4
    iterOpts = struct('Debug',false, 'NumCores', 4);
else
    p              = inputParser();
    p.StructExpand = true;
    addParameter(p,'Debug',false);
    addParameter(p,'NumCores',4);
    p.parse(iterOpts);
    iterOpts       = p.Results;
end

% Pass median prices
iterOpts.CacheByFile = loadlatest('median');
iterOpts.CacheByFile = cache2cell(iterOpts.CacheByFile, iterOpts.CacheByFile.File);

% Single file or whole datastore
if nargin == 1 && isstruct(path2data)
    res = getBadPricesCountSingle_(path2data, [], [], cleanOpts);
else
    res = getBadPricesCountWholeDatastore_(path2data, outname, cleanOpts, iterOpts);
end
end

function res = getBadPricesCountSingle_(s, cached, filenum, cleanOpts)
[~, iexclude, nobs] = clean_consolidate(s, cached, cleanOpts);

% Counts
res         = s.index(:,{'Id','Date','Permno'});
subs        = uint32(RunLength((1:size(cached,1))',nobs));
res.Ntot    = uint32(nobs);
res.Nbadtot = uint32(accumarray(subs,  iexclude));
end

function res = getBadPricesCountWholeDatastore_(path2data, outname, cleanOpts, iterOpts)
iterOpts = [fieldnames(iterOpts), struct2cell(iterOpts)]';
res      = taq.iterate_datastore(path2data, outname,...
    @(s,cached,filenum) getBadPricesCountSingle_(s, cached, filenum, cleanOpts), iterOpts{:});
end
