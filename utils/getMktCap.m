function cap = getMktCap(permno, date, iscommon,issp500, ispanel, lag)
if nargin < 1 || isempty(permno),   permno   = [];      end
if nargin < 2 || isempty(date),     date     = [];      end
if nargin < 3 || isempty(iscommon), iscommon = true;    end
if nargin < 4 || isempty(issp500),  issp500  = false;   end
if nargin < 5 || isempty(ispanel),  ispanel  = false;   end

try
    cap = loadresults('mktcap');
catch
    cap = loadresults('mktcap','..\results');
end

if nargin == 6 && lag > 0
    % Ensure it is sorted by id-date, i.e. some date changes correspond 
    % to same permno, we need to sort
    idx = diff(cap.Date) == 0;
    if ~all(cap.Permno(idx) ~= cap.Permno([false;idx]))
        cap = sortrows(cap,{'Permno','Date'});
    end
    
    % Lag
    cap.Cap = [NaN(lag,1); cap.Cap(1:end-lag)];
    idx     = [false(lag,1); cap.Permno(1:end-lag) == cap.Permno(1+lag:end)];
    cap     = cap(idx,:);
end

% Filter by permno, date
idx = true(size(cap,1),1);
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