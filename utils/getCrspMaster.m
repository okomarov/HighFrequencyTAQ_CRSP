function record = getCrspMaster(idtype,id,date)

master = loadresults('msenames');

% CRSP
tmp = master(ismember(master.(upperfirst(idtype)), id),:);
tmp = sortrows(tmp, 'Namedt');
if size(tmp,1) == 1
    record = tmp;
else
    record = tmp(in(date,[tmp.Namedt, tmp.Nameendt]),:);
end
record.Date = date;
end