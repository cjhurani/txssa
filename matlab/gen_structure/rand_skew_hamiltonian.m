function H = rand_skew_hamiltonian(n)

% H = rand_skew_hamiltonian(n)
%
% Output a random dense real  skew-Hamiltonian matrix H of size n.  n must be
% even.

assert(mod(n,2) == 0, 'rand_skew_hamiltonian: input n = %d must be even', n);

n = n/2;

B = zrand(n) - 0.5;
C = zrand(n) - 0.5; C = 0.5*(C - C');
D = zrand(n) - 0.5; D = 0.5*(D - D');

H = [B C; D B'];

