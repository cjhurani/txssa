function A = rand_symmetric(n)

% A = rand_symmetric(n)
%
% Output a random dense real symmetric matrix A of size n.

A = zrand(n) - 0.5;

A = 0.5*(A + A');


