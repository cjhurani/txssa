function C = rand_skew_circulant_complex(n)

% C = rand_skew_circulant_complex(n)
%
% Output a random dense complex  skew_circulant matrix C of size n.

C = skew_circulant(complex(zrand(n,1) - 0.5, zrand(n,1) - 0.5));

