function res = estimateMktCap()
% Load data
id2permno = loadresults('masterPermno');

dsfquery     = loadresults('dsfquery');
dsfquery     = dsfquery(:,{'Permno', 'Date','Prc'});
dsfquery.Prc = abs(dsfquery.Prc);

dseshares = loadresults('dseshares');
dseshares = unique(dseshares(:,{'PERMNO', 'SHRSDT','SHRENDDT','SHROUT'}));

% Reference permnos and dates
Permno = intersect(dsfquery.Permno, dseshares.PERMNO);
Permno = intersect(Permno, id2permno.Permno);
Date   = unique(id2permno.Date);
sz     = [numel(Date), numel(Permno)];

% Filter
idx      = ismember(dsfquery.Permno,Permno) & ismember(dsfquery.Date,Date);
dsfquery = dsfquery(idx,:);

idx       = ismember(dseshares.PERMNO, Permno) & dseshares.SHRENDDT >= min(Date);
dseshares = dseshares(idx,:);

% Unstack prices
[~,row]        = ismember(dsfquery.Date, Date);
[~,col]        = ismember(dsfquery.Permno, Permno);
dsfquery       = accumarray([row,col], dsfquery.Prc, sz);
inan           = ~accumarray([row,col], 1, sz);
dsfquery(inan) = NaN;
dsfquery       = nanfillts(dsfquery,true);

% Unstack number of shares
dseshares        = convertColumn(dseshares,'double','SHROUT');
dseshares        = pivotFromTo(dseshares);
irow             = ismember(dseshares.Date, Date);
icol             = [false; ismember(dseshares.Id, Permno)];
dseshares        = dseshares.Panel(irow,icol);
dseshares        = dseshares{:,:};
inull            = dseshares == 0;
dseshares(inull) = NaN;
dseshares        = nanfillts(dseshares);

% Market capitalizations
res = dseshares.*dsfquery;

% Unstack
[Date,Permno] = ndgrid(Date,Permno);
notnan        = ~isnan(res);
res           = table(Date(notnan), Permno(notnan), res(notnan),...
                'VariableNames',{'Date','Permno','Cap'});
            
filename = sprintf('%s_mktcap.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile('results',filename), 'res')
end