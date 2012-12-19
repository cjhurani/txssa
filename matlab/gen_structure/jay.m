function J = jay(n)

% J = jay(n)
%
% Creates a dense matrix J with ones on the antidiagonal.

J = fliplr(eye(n));

