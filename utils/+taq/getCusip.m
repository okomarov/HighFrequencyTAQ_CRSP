function cusip = getCusip(symbol,date,path2mst)
% GETCUSIP Retrieves CUSIPs of a symbol
%
%   GETCUSIP(SYMBOL, DATE)
%       Retrieves the CUSIPs corresponding to the SYMBOL - DATE(s) pairs.
%
% See also: getMaster
if ~isrowchar(symbol)
    error('taq:getCusip:invalidSymbol','SYMBOL should be a string.')
end
if nargin < 3 || isempty(path2mst)
    path2mst = '.\data\TAQ\master';
end
dttype = getDateType(date);

switch dttype
    case {'scalar','set'}
        mst = taq.getMaster(symbol,path2mst);
        mst = mst{1};

        if size(mst.CUSIP,2) < 8
            cusip = '';
            return
        end

        [~,isort] = sort(mst.FDATE);
        mst       = mst(isort,:);

        [~,~,bin] = histcounts(date, [mst.FDATE; inf]);
        cusip     = mst.CUSIP(bin,1:8);

    otherwise
        error('taq:getCusip:invalidDateType','Only scalar or a set of yyyymmdd dates are accepted.')
end
end
