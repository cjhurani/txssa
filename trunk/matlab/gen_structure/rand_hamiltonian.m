function H = rand_hamiltonian(n)

% H = rand_hamiltonian(n)
%
% Output a random dense real  Hamiltonian matrix H of size n.  n must be even.

assert(mod(n,2) == 0, 'rand_hamiltonian: input n = %d must be even', n);

n = n/2;

A = zrand(n) - 0.5;
B = zrand(n) - 0.5; B = 0.5*(B + B');
C = zrand(n) - 0.5; C = 0.5*(C + C');

H = [A B; C (-A')];
