function tb = getMktCap(tb, lag, ispanel)
if nargin < 2,                      lag     = [];      end
if nargin < 3 || isempty(ispanel),  ispanel = false;   end

try
    cap = loadresults('mktcap');
catch
    cap = loadresults('mktcap','..\results');
end

try
    cap = lagpanel(cap,'Permno',lag);
catch
    cap = sortrows(cap,{'Permno','Date'});
    cap = lagpanel(cap,'Permno',lag);
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
