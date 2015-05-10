function tb = consolidateFromTo(tb)
% CONSOLIDATEFROMTO Remove intermediate redundant from-to records
%
%   CONSOLIDATEFROMTO(TB) TB should be a table with:
%           ID | From | To | Value

% Track changes in the value, i.e. Id | VALUE | FROM 
idx = isfeatchange(tb(:,[1,4,2]));
% Start/end positions
stpos   = find(idx);
enpos   = [stpos(2:end)-1; size(tb,1)];
% Enddates
enddate = tb(enpos,3); 
% Filter out
tb      = tb(stpos,:);
% Assign enddates back
tb(:,3) = enddate;
end