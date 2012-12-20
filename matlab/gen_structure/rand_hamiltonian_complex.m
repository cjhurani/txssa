function H = rand_hamiltonian_complex(n)

% H = rand_hamiltonian_complex(n)
%
% Output a random dense complex  Hamiltonian matrix H of size n.  n must be even.

assert(mod(n,2) == 0, 'rand_hamiltonian_complex: input n = %d must be even', n);

n = n/2;

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);
B = complex(zrand(n) - 0.5, zrand(n) - 0.5); B = 0.5*(B + B');
C = complex(zrand(n) - 0.5, zrand(n) - 0.5); C = 0.5*(C + C');

H = [A B; C (-A')];
