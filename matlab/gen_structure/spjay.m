function J = spjay(n)

% J = spjay(n)
%
% Creates a sparse matrix J with ones on the antidiagonal.

J = fliplr(speye(n));

