function H = rand_skew_hamiltonian_complex(n)

% H = rand_skew_hamiltonian_complex(n)
%
% Output a random dense complex  skew-Hamiltonian matrix H of size n.  n must be
% even.

assert(mod(n,2) == 0, 'rand_skew_hamiltonian_complex: input n = %d must be even', n);

n = n/2;

B = complex(zrand(n) - 0.5, zrand(n) - 0.5);
C = complex(zrand(n) - 0.5, zrand(n) - 0.5); C = 0.5*(C - C');
D = complex(zrand(n) - 0.5, zrand(n) - 0.5); D = 0.5*(D - D');

H = [B C; D B'];

