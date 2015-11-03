function cached = cache2cell(data, groupvar, sortMode)
% [cached,unGroups] = cache2cell(data, groupvar, sortmode)
%
% NOTE: it sorts data according to groupvar and operates along rows
if size(data,1) ~= numel(groupvar)
    error('cache2cell:sizeMismatch','DATA must have as many rows as elements in the GROUPVAR.')
end
if nargin == 3
    [groupvar,isort] = sort(groupvar,[],sortMode);
    data             = data(isort,:);
end
cached = accumarray(groupvar,(1:size(data))',[],@(x) {data(x,:)});
end