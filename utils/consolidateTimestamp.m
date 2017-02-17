function data = consolidateTimestamp(data, type)
% CONSOLIDATETIMESTAMP Consolidates prices and volume with same timestamp
%
%   CONSOLIDATETIMESTAMP(DATA, [TYPE])
%       - DATA is a table with Id or Permno, Datetime or Date and Time,
%         Price and (optionally) Volume fields
%       - TYPE can be 'volumeWeighted' (default), 'first', 'median' or 'skip'
%
%   NOTE: mostly used on TAQ trades data

if nargin < 2
    type = 'volumeWeighted'; 
end

try
    dt = data.Datetime;
catch
    dt = hhmmssmat2serial(data.Time) + yyyymmdd2serial(data.Date);
end
try
    id = double(data.Id);
catch
    id = double(data.Permno);
end
[~,pos,subs] = unique([id,dt],'rows');
try
    vol    = double(data.Volume)/100;
    voltot = accumarray(subs, vol);
catch
    type = 'median';
    warning('DATA has no ''Volume'' field. Using ''median'' consolidation.')
end

switch type
    case 'first'
        idx    = [true; logical(diff(subs))];
        prices = data.Price(idx);
    case 'median'
        prices = accumarray(subs, data.Price,[],@fast_median);
    case 'volumeWeighted'
        prices = accumarray(subs, data.Price.*vol) ./ voltot;
    case 'skip'
        % Do not consolidate price
end
data       = data(pos,:);
data.Price = prices;
try
    data.Volume = voltot*100;
catch
end
end