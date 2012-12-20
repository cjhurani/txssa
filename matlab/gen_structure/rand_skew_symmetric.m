function A = rand_skew_symmetric(n)

% A = rand_skew_symmetric(n)
%
% Output a random dense real skew_symmetric matrix A of size n.

A = zrand(n) - 0.5;

A = 0.5*(A - A');

