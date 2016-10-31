function out = getCusip(symbol,date,path2mst)
% GETCUSIP Retrieves CUSIPs of a symbol
%
%   GETCUSIP(SYMBOL, [DATE])
%       Retrieves all CUSIPs corresponding to a SYMBOL or a
%       cell array of symbols. Can filter by DATE (default: takes all).
%
% See also: getMaster


if nargin < 2 || isempty(date),        date     = inf;                 end
if nargin < 3 || isempty(path2mst),    path2mst = '.\data\TAQ\master'; end

allMst = taq.getMaster(symbol,path2mst);
dttype = getDateType(date);

n   = numel(symbol);
out = cell(n,1);

for ii = 1:n
    mst = allMst{ii};
    if size(mst.CUSIP,2) < 8
        continue
    end

    mst.CUSIP = mst.CUSIP(:,1:8);
    mst = sortrows(mst(:,{'SYMBOL','CUSIP','FDATE'}),'FDATE');

    [idx,bin] = filterByDate(mst.FDATE, date, dttype);

    if strcmpi(dttype, 'set')
        out{ii} = table(mst.SYMBOL(bin,:), mst.CUSIP(bin,:), date,...
            'VariableNames',{'SYMBOL','CUSIP','FDATE'});
    else
        % Reduce records to earliest
        mst     = mst(idx,:);
        ikeep   = isfeatchange(mst(:,[2,3]), 2);
        out{ii} = mst(ikeep,:);
    end
end
end

function [idx,bin] = filterByDate(dates, dt, type)
nrows = size(dates,1);
switch type
    case {'scalar','set'}
        [~,~,bin] = histcounts(dt, [dates; inf]);
        idx       = false(nrows,1);
        idx(bin)  = true;
        
    case 'all'
        idx = true(nrows,1);
        
    case 'ge'
        idx = dates > dt(1);
        pos = find(idx,1,'first');
        if pos > 1
            idx(pos-1) = true;
        end
        
    case 'le'
        idx = dates <= dt(2);
        
    case 'fromto'
        idx = dates >= dt(1) & dates <= dt(2);

end
end
