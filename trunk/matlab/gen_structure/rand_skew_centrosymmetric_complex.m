function R = rand_skew_centrosymmetric_complex(n)

% R = rand_skew_centrosymmetric_complex(n)
%
% Output a random dense complex  skew-centrosymmetric matrix R of size n.

if(mod(n,2) == 0)
    n = n/2;
    even = true;
else
    n = (n-1)/2;
    even = false;
end

A = complex(zrand(n) - 0.5, zrand(n) - 0.5);
B = complex(zrand(n) - 0.5, zrand(n) - 0.5);
J = spjay(n);

if(even)
    R = [A (-B*J); J*B (-J*A*J)];
else
    x = complex(zrand(n,1) - 0.5, zrand(n,1) - 0.5);
    y = complex(zrand(1,n) - 0.5, zrand(1,n) - 0.5);
    R = [A J*x (-B*J); y*J 0 (-y); J*B (-x) (-J*A*J)];
end

