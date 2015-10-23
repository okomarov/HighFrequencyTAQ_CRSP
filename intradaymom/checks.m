%% Check Open/Close in TAQ vs CRSP
taq = loadresults('price_fl');

if OPT_NOMICRO
    idx = isMicrocap(taq,'LastPrice',OPT_LAGDAY);
    taq = taq(~idx,:);
end

crsp      = loadresults('dsfquery');
crsp.Prc  = abs(crsp.Prc);
[~,ia,ib] = intersectIdDate(crsp.Permno,crsp.Date, taq.Permno, taq.Date);
crsp      = crsp(ia,:);
taq       = taq(ib,:);

% Filter out outliers
taq.TAQret   = taq.LastPrice./taq.FirstPrice-1;
iout         = taq.TAQret          > OPT_OUTLIERS_THRESHOLD |...
               1./(taq.TAQret+1)-1 > OPT_OUTLIERS_THRESHOLD;
taq(iout,:)  = [];            
crsp(iout,:) = [];
            
% Comparison table
cmp         = [crsp(:,{'Date','Permno','Openprc'}) taq(:,'FirstPrice'), ...
    crsp(:,{'Bid','Ask','Prc'}),taq(:,{'LastPrice','TAQret'})];
cmp.CRSPret = cmp.Prc./cmp.Openprc-1;

retdiff = abs(cmp.CRSPret - cmp.TAQret);
idx     = retdiff > eps*1e12;
boxplot(retdiff(idx))
