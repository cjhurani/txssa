function K = kay(n)

% K = kay(n)
%
% Creates a dense K = [0 I; -I 0] of size n if n is even, else error.

assert(mod(n,2) == 0, 'kay: input n = %d must be even', n);

n = n/2;
z = zeros(n);
i = eye(n);
K = [z i; -i z];
