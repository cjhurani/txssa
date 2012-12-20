function R = rand_skew_persymmetric(n)

% R = rand_skew_persymmetric(n)
%
% Output a random dense real  skew_persymmetric matrix R of size n.

A = zrand(n) - 0.5;
J = spjay(n);

R = 0.5*(A - J*A'*J);
