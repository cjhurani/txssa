function A = rand_skew_complex_symmetric(n)

% A = rand_skew_complex_symmetric(n)
%
% Output a random dense skew_complex_symmetric matrix A of size n.

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);

A = 0.5*(A - A.');

