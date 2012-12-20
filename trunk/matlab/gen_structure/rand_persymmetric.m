function R = rand_persymmetric(n)

% R = rand_persymmetric(n)
%
% Output a random dense real  persymmetric matrix R of size n.

A = zrand(n) - 0.5;
J = spjay(n);

R = 0.5*(A + J*A'*J);
