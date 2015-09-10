function A = nonnan(A)
A = A(~isnan(A));
end