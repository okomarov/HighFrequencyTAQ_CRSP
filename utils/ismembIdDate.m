function [idx, pos, keyA, keyB] = ismembIdDate(IdA, DateA, IdB, DateB)
if numel(IdA) ~= numel(DateA)
    error('ismembIdDate:numElem','ID and DATE for A must have same length.')
end
if numel(IdB) ~= numel(DateB)
    error('ismembIdDate:numElem','ID and DATE for B must have same length.')
end
% Build composite key from id and date
keyA       = uint64(IdA)*1e8 + uint64(DateA);
keyB       = uint64(IdB)*1e8 + uint64(DateB);
[idx, pos] = ismember(keyA,keyB);
end