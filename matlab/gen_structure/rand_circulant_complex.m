function C = rand_circulant_complex(n)

% C = rand_circulant_complex(n)
%
% Output a random dense complex  circulant matrix C of size n.

C = circulant(complex(zrand(n,1) - 0.5, zrand(n,1) - 0.5));

