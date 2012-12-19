function C = rand_circulant(n)

% C = rand_circulant(n)
%
% Output a random dense real  circulant matrix C of size n.

C = circulant(zrand(n,1) - 0.5);

