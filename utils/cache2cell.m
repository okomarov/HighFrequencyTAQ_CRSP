function [cached,unGroups] = cache2cell(data, groupvar)
% [cached,unGroups] = cache2cell(data, groupvar)
%
% NOTE: it sorts data according to groupvar and operates along rows

if size(data,1) ~= numel(groupvar)
    error('cache2cell:sizeMismatch','DATA must have as many rows as elements in the GROUPVAR.')
end
if ~issorted(groupvar)
    [groupvar,isort]  = sort(groupvar);
    data              = data(isort,:);
end
[unGroups,~,subs] = unique(groupvar);
nrows             = accumarray(subs,1);
cached            = mat2cell(data,nrows,size(data,2));
end