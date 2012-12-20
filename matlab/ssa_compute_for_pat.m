function X = ssa_compute_for_pat(A, A_pat, max_num_bins, impose_null_spaces)

% X = ssa_compute(A, ratio, p, max_num_bins, impose_null_spaces)
%
% Wrapper function over all the sparse spectral approximation (ssa)
% functionality.
%
% Given
%
%   a matrix A
%   a pattern A_pat of same size as A with zeros and ones
%   a non-negative integer max_num_bins (200-1000 is a reasonable choice)
%   an option to impose null-spaces     (true or false)
%
% this function computes a matrix X that is sparse (like A_pat is)
% and spectrally close to A, especially in the lower end of the singular
% value spectrum.  Even if A is complex, A_pat must be real.
%
% If impose_null_spaces is false, returned without imposing null space.  If
% impose_null_spaces is true, null-spaces (if any) are imposed on X.
% max_num_bins == 0 means no binning is performed and can lead to significant
% slowdown.
%
% No error checking is performed in this function.  It is assumed that
% the functions that are called will do the checks.

[pinv_A, rnull, lnull] = pinv_rrqr(A);

A_id = bin_sparse_matrix(A, A_pat, max_num_bins);

pinv_ATA = pinv_A * pinv_A';
pinv_AAT = pinv_A' * pinv_A;

X = ssa_minimization(A, A_id, pinv_ATA, pinv_AAT);

if(impose_null_spaces && (0 < size(rnull, 2) || 0 < size(lnull, 2)))
    X = ssa_impose_action(X, A, rnull, lnull);
end
