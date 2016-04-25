function permno = getPermno(ncusip)
% GETPERMNO Retrieve CRSP Permno from historical CUSIP
%
% See also: getNames

names       = crsp.getNames('monthly');
[idx,pos]   = ismember(cellstr(ncusip), names.Ncusip);
permno      = zeros(size(ncusip,1),1,'uint32');
permno(idx) = names.Permno(pos(idx));
end
