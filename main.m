% Num issues matched
path2data = '.\data\TAQ';
master    = load(fullfile(path2data, 'master'), '-mat');
master    = addPermno(master.mst);

[dates,~,subs] = unique(master.Date);
nobs           = master.To-master.From+1;
tot            = accumarray(subs,nobs);
idx            = master.Permno~=0;
matched        = accumarray(subs(idx),nobs(idx));

mean(matched./tot)
min(matched./tot)