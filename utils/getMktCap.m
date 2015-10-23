function tb = getMktCap(tb, lag, ispanel)
if nargin < 3 || isempty(ispanel),  ispanel = false;   end

try
    cap = loadresults('mktcap');
catch
    cap = loadresults('mktcap','..\results');
end

if nargin >= 2 && lag > 0
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
if isnumeric(tb)
    tb = cap(ismember(cap.Permno, tb),:);
else
    [idx,pos]     = ismembIdDate(tb.Permno, tb.Date, cap.Permno, cap.Date);
    tb.Cap(idx,1) = cap.Cap(pos(idx));
end

% Unstack
if ispanel
    tb = sortrows(unstack(tb(:,{'Cap','Permno','Date'}),'Cap','Permno'),'Date');
end
end