function J = ssa_misfit(A, X)

pA = pinv_rrqr(A);

J = 0.5 * (norm((X - A) * pA, 'fro')^2 + norm(pA * (X - A), 'fro')^2);

