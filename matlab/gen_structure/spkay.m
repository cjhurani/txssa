function K = spkay(n)

% K = spkay(n)
%
% Creates a sparse K = [0 I; -I 0] of size n if n is even, else error.

assert(mod(n,2) == 0, 'spkay: input n = %d must be even', n);

n = n/2;
z = sparse(n,n);
i = speye(n);
K = [z i; -i z];
