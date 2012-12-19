function pc = pcond(A)

s = sort(svd(full(A)));
tol = 100 * min(size(A)) * eps(max(s));  % MAGIC CONSTANT
r = sum(s > tol);
pc = s(end)/s(size(s,1) - r + 1);
