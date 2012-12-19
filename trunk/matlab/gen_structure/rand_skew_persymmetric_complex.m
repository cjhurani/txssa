function R = rand_skew_persymmetric_complex(n)

% R = rand_skew_persymmetric_complex(n)
%
% Output a random dense complex  skew_persymmetric matrix R of size n.

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);
J = spjay(n);

R = 0.5*(A - J*A'*J);
