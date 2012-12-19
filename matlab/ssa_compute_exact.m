function X = ssa_compute_exact(A, ratio, p)

% X = ssa_compute_exact(A, ratio, p)
%
% Wrapper function over all the sparse spectral approximation (ssa)
% functionality.  This one uses the "exact" algorithm which means
% no binning is done, and null-space constraints are imposed while
% minimization rather than in the last step.
%
% Given
%
%   a matrix A
%   a sparsity ratio in [0,1]           (0 means more sparse, 1 means less)
%   a power p in [0,Inf]                (for L_p norm based sparsity.  1 is OK.)
%
% this function computes a matrix X that is sparse and spectrally close to A,
% especially in the lower end of the singular value spectrum.
%
% No error checking is performed in this function.  It is assumed that
% the functions that are called will do the checks.

[pinv_A, rnull, lnull] = pinv_rrqr(A);

%fringe case: If a row is zero, the left null-space will
%have a "trivial" component.  They should not be needed in
%min_per_row.  Same for columns.

[num_near_zero_rows num_near_zero_cols] = near_zero_row_col(A);

min_per_row = max(0,min(size(rnull, 2) - num_near_zero_cols, size(A,2)));
min_per_col = max(0,min(size(lnull, 2) - num_near_zero_rows, size(A,1)));
A_pat = p_norm_sparsity_matrix(A, ratio, p, min_per_row, min_per_col);

pinv_ATA = pinv_A * pinv_A';
pinv_AAT = pinv_A' * pinv_A;

X = ssa_minimization_with_action(A, A_pat, pinv_ATA, pinv_AAT, rnull, lnull);
