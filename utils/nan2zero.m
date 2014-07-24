function A = nan2zero(A)
% NAN2ZERO Replaces NaNs with 0s
%
%   NAN2ZERO(A)
A(isnan(A)) = 0;
end