function C = circulant(c)

% C = circulant(c)
%
% Output a circulant matrix C with its first column equal to the vector
% given by c.  The input can be real or complex.

assert(isvector(c), 'circulant: input array must be a 1-D vector');

m = length(c);
c = c(:);
idx = [m 1:(m-1)];

C = zeros(m);

if(0 < m)
    C(:,1) = c;
end

for j = 2:m
    C(:,j) = C(idx,j-1);
end
