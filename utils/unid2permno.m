function permnos = unid2permno(unids)
% UNID2PERMNO Returns the permno corresponding to the UnID
%
%   UNID2PERMNO(UNIDS) Array of UnIDs
%
%
%   PERMNOS = UNID2PERMNO(...) Array of PERMNOS in preserved order.
%
% Note: only one permno per UnID (excluding NaN permnos).


taq2crsp = loadresults('taq2crsp');
taq2crsp = unique(taq2crsp(:,{'permno','ID'}));
[~,pos]  = ismember(unids, taq2crsp.ID);
permnos  = taq2crsp.permno(pos);

end