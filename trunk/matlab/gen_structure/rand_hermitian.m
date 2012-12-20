function A = rand_hermitian(n)

% A = rand_hermitian(n)
%
% Output a random dense complex hermitian matrix A of size n.

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);

A = 0.5*(A + A');


