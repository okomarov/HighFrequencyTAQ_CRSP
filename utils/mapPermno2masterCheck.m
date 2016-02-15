mstp   = loadresults('masterPermno');
master = load(fullfile('.\data\TAQ','master'),'-mat','ids');
ids    = master.ids;
master = loadresults('TAQmaster');
mse    = loadresults('msenames');

% TAQ
permno = 0;
while permno == 0
    pos    = randi(size(mstp,1));
    permno = mstp.Permno(pos);
end
date   = mstp.Date(pos);
id     = mstp.Id(pos);
tmp    = master(ismember(master.SYMBOL,ids{id}),:);
tmp    = sortrows(tmp, 'FDATE');
if size(tmp,1) == 1
    taq = tmp(:,1:3);
else
    [~,~,bin] = histcounts(date,[double(tmp.FDATE)-0.5; inf]);
    taq    = tmp(bin,1:3);
end
taq.Date = date;
taq.Id   = id;
taq.Permno = permno;

% CRSP
crsp = mse(mse.Permno == taq.Permno,:);
crsp = crsp(in(date, [crsp.Namedt, crsp.Nameendt]),:);
crsp = crsp(:,[11,9,7,2,3,1]);