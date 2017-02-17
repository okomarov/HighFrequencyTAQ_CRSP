function res = getMedian(path2data, outname, cleanOpts, iterOpts)
% Parse options
if nargin < 3
    cleanOpts = struct('ExcludeBadDays',false, 'BadPriceMultiplier', [],'ConsolidateTimestampType',[]);
else
    p              = inputParser();
    p.StructExpand = true;
    addParameter(p,'BadPriceMultiplier',[]);
    addParameter(p,'ExcludeBadDays',false);
    addParameter(p,'ConsolidateTimestampType',[]);
    p.parse(cleanOpts);
    cleanOpts      = p.Results;
end

% Single file or whole datastore
if nargin == 1 && isstruct(path2data)
    res = getMedianSingle_(path2data, [], [], cleanOpts);
else
    res = getMedianWholeDatastore_(path2data, outname, cleanOpts);
end
end

function res = getMedianSingle_(s, cached, filenum, cleanOpts)
data = clean_consolidate(s, cached, cleanOpts);

[un,~,subs] = unique([uint32(data.Id),data.Date],'rows');
nmst        = size(un,1);
prices      = accumarray(subs, data.Price,[nmst,1], @fast_median);

res                    = s.index(:,{'Id','Date','Permno'});
res.MedianPrice(:,1)   = NaN;
[~,pos]                = ismembIdDate(un(:,1), un(:,2), res.Id, res.Date);
res.MedianPrice(pos,1) = prices;
end

function res = getMedianWholeDatastore_(path2data, outname, cleanOpts)
res = taq.statsByDay(path2data, outname, @(s,cached,filenum) getMedianSingle_(s, cached, filenum, cleanOpts));
end
