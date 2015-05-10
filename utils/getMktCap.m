function cap = getMktCap(permno, date, iscommon,issp500, ispanel)
if nargin < 1 || isempty(permno),   permno   = [];      end
if nargin < 2 || isempty(date),     date     = [];      end
if nargin < 3 || isempty(iscommon), iscommon = true;    end
if nargin < 4 || isempty(issp500),  issp500  = false;   end
if nargin < 5 || isempty(ispanel),  ispanel  = false;   end

cap = loadresults('mktcap');
idx = true(size(cap,1),1);

% Filter by permno, date
if ~isempty(permno) && ~isempty(date)
    idx = idx & ismembIdDate(cap.Permno,cap.Date,permno,date);
else
    if ~isempty(permno)
        idx = idx & ismember(cap.Permno,permno);
    end
    if ~isempty(date)
        idx = idx & ismember(cap.Date,date);
    end
end

% Common share and sp500 constituents
if iscommon
    idx = idx & iscommonshare(cap);
end
if issp500
    idx = idx & issp500member(cap);
end
cap = cap(idx,:);

% Unstack
if ispanel
    cap = sortrows(unstack(cap,'Cap','Permno'),'Date');
end
end