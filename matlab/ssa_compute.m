function X = ssa_compute(A, ratio, p, max_num_bins, impose_null_spaces)

% X = ssa_compute(A, ratio, p, max_num_bins, impose_null_spaces)
%
% Wrapper function over all the sparse spectral approximation (ssa)
% functionality.
%
% Given
%
%   a matrix A
%   a sparsity ratio in [0,1]           (0 means more sparse, 1 means less)
%   a power p in [0,Inf]                (for L_p norm based sparsity.  1 is OK.)
%   a non-negative integer max_num_bins (200-1000 is a reasonable choice)
%   an option to impose null-spaces     (true or false)
%
% this function computes a matrix X that is sparse and spectrally close to A,
% especially in the lower end of the singular value spectrum.
%
% If impose_null_spaces is false, returned without imposing null space.  If
% impose_null_spaces is true, null-spaces (if any) are imposed on X.
% max_num_bins == 0 means no binning is performed and can lead to significant
% slowdown.
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

A_id = bin_sparse_matrix(A, A_pat, max_num_bins);

pinv_ATA = pinv_A * pinv_A';
pinv_AAT = pinv_A' * pinv_A;

X = ssa_minimization(A, A_id, pinv_ATA, pinv_AAT);

if(impose_null_spaces && (0 < size(rnull, 2) || 0 < size(lnull, 2)))
    X = ssa_impose_action(X, A, rnull, lnull);
end
