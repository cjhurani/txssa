function A = rand_complex_symmetric(n)

% A = rand_complex_symmetric(n)
%
% Output a random dense complex_symmetric matrix R of size n.

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);

A = 0.5*(A + A.');

