function R = rand_skew_centrosymmetric(n)

% R = rand_skew_centrosymmetric(n)
%
% Output a random dense real  skew-centrosymmetric matrix R of size n.

if(mod(n,2) == 0)
    n = n/2;
    even = true;
else
    n = (n-1)/2;
    even = false;
end

A = zrand(n) - 0.5;
B = zrand(n) - 0.5;
J = spjay(n);

if(even)
    R = [A (-B*J); J*B (-J*A*J)];
else
    x = zrand(n,1) - 0.5;
    y = zrand(1,n) - 0.5;
    R = [A J*x (-B*J); y*J 0 (-y); J*B (-x) (-J*A*J)];
end
