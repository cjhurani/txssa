function C = rand_skew_circulant(n)

% C = rand_skew_circulant(n)
%
% Output a random dense real  skew_circulant matrix C of size n.

C = skew_circulant(zrand(n,1) - 0.5);

