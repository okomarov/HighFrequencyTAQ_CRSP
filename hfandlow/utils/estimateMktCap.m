function res = estimateMktCap()

from = 19930101;

% Load prices
dsfquery     = loadresults('dsfquery');
dsfquery     = dsfquery(:,{'Permno', 'Date','Prc'});
dsfquery.Prc = abs(dsfquery.Prc);
idx          = dsfquery.Date >= from;
dsfquery     = dsfquery(idx,:);
to           = max(dsfquery.Date);

% Load shares
dseshares               = loadresults('dseshares');
dseshares               = unique(dseshares(:,{'PERMNO', 'SHRSDT','SHRENDDT','SHROUT'}));
idx                     = dseshares.SHRENDDT >= from;
dseshares               = dseshares(idx,:);
idx                     = dseshares.SHRENDDT > to;
dseshares.SHRENDDT(idx) = to;
idx                     = dseshares.SHRSDT < from;
dseshares.SHRSDT(idx)   = from;

% Intersect permnos
Permno    = intersect(dsfquery.Permno, dseshares.PERMNO);
idx       = ismember(dsfquery.Permno,Permno);
dsfquery  = dsfquery(idx,:);
idx       = ismember(dseshares.PERMNO, Permno);
dseshares = dseshares(idx,:);

Date = unique(dsfquery.Date);

% Unstack prices
dsfquery = sortrows(tbextend.unstack(dsfquery,'Prc','Permno'),'Date');
dsfquery = nanfillts(dsfquery{:,2:end},true);

% Unstack number of shares
dseshares        = consolidateFromTo(dseshares);
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

% Stack back
[Date,Permno] = ndgrid(Date,Permno);
notnan        = ~isnan(res);
res           = table(Date(notnan), Permno(notnan), res(notnan),...
                'VariableNames',{'Date','Permno','Cap'});
            
filename = sprintf('%s_mktcap.mat',datestr(now,'yyyymmdd_HHMM'));
save(fullfile('results',filename), 'res')
end