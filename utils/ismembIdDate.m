function [idx, pos, keyA, keyB] = ismembIdDate(IdA, DateA, IdB, DateB)
% Build composite key from id and date
keyA       = uint64(IdA)*1e8 + uint64(DateA);
keyB       = uint64(IdB)*1e8 + uint64(DateB);
[idx, pos] = ismember(keyA,keyB);
end