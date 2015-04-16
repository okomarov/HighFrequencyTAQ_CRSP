function idx = isfeatchange(tb, inum)

% ISFEATCHANGE Progressively detects changes across features/ids
%
%   ... = ISFEATCHANGE(TB) TB should be a table with: 
% 
%               ID | feature cols | DATES
%          
%         ISFEATCHANGE(..., INUM) INUM is a logical vector telling
%                                 which columns are numeric. By default
%                                 INUM is all true. 
%                                 Non numeric columns are considered 
%                                 cellstrings.
%
%   IDX = ... 
%         Returns the index to the records to keep
%
%
%   NOTE: this function is used to retain the earliest-in-time records when
%         a feture changed. It is designed to avoid accumarray() in the
%         following situation:
%           
%         [Input]              [accumarray(...,@min(date))]            
%
%         feat         date    
%            2   2012/12/22     1
%            2   2012/12/31     0
%            1   2013/01/05     1
%            2   2013/02/07     0 <- feature change not picked!

% Properly sorted
[~,~,subs] = unique(tb(:,[1,end]));
if ~issorted(subs)
    error('tconsolidate:unsorted', 'TB should be sorted by Id and Date.')
end
% Inum
szTb = size(tb);
if nargin < 2 || isempty(inum)
    inum = 1:szTb(2);
elseif islogical(inum) 
    inum = find(inum);
end

vnames = getVariableNames(tb);
idx = false(size(tb,1)-1,1);
for c = 1:szTb(2)-1
    field = vnames{c};
    if any(c == inum)
        idx = idx | tb.(field)(2:end) ~= tb.(field)(1:end-1);
    else
        idx = idx | ~strcmpi(tb.(field)(2:end),tb.(field)(1:end-1));
    end
end
idx = [true; idx];
end