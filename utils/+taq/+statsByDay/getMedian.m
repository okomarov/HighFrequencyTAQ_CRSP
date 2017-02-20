function res = getMedian(path2data, outname, cleanOpts, iterOpts)
% Parse options
if nargin < 3 || isempty(cleanOpts)
    cleanOpts = struct('ExcludeBadDays',false, 'BadPriceMultiplier', [], 'TerminateEarly',false,'ConsolidateTimestampType',[]);
else
    p              = inputParser();
    p.StructExpand = true;
    addParameter(p,'BadPriceMultiplier',[]);
    addParameter(p,'ExcludeBadDays',false);
    addParameter(p,'TerminateEarly',false);
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

% Single file or whole datastore
if nargin == 1 && isstruct(path2data)
    res = getMedianSingle_(path2data, [], [], cleanOpts);
else
    res = getMedianWholeDatastore_(path2data, outname, cleanOpts, iterOpts);
end
end

function res = getMedianSingle_(s, cached, filenum, cleanOpts)
data = clean_consolidate(s, cached, cleanOpts);

[un,~,subs] = unique([uint32(data.Id),data.Date],'rows');
nmst        = size(un,1);
prices      = single(accumarray(subs, data.Price,[nmst,1], @fast_median));

res                    = s.index(:,{'Id','Date','Permno'});
res.MedianPrice(:,1)   = NaN;
[~,pos]                = ismembIdDate(un(:,1), un(:,2), res.Id, res.Date);
res.MedianPrice(pos,1) = prices;
res.File(:,1)          = uint16(filenum);
end

function res = getMedianWholeDatastore_(path2data, outname, cleanOpts, iterOpts)
iterOpts = [fieldnames(iterOpts), struct2cell(iterOpts)]';
res      = taq.iterate_datastore(path2data, outname, @(s,cached,filenum) getMedianSingle_(s, cached, filenum, cleanOpts), iterOpts{:});
end
