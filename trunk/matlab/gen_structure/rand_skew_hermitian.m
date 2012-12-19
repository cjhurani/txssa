function A = rand_skew_hermitian(n)

% A = rand_skew_hermitian(n)
%
% Output a random dense complex skew_hermitian matrix A of size n.

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);

A = 0.5*(A - A');



