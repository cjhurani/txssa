function C = skew_circulant(c)

% C = skew_circulant(c)
%
% Output a skew_circulant matrix C with its first column equal to the vector
% given by c.  The input can be real or complex.

assert(isvector(c), 'skew_circulant: input array must be a 1-D vector');

C = circulant(c);
C = C .* (1 - 2 * triu(ones(length(c)),1));
