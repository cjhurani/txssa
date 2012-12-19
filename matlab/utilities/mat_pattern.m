function pat = mat_pattern(A)

pat = A;
pat(pat ~= 0) = 1;
pat = sparse(pat);
