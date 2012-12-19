function X = ssa_compute_exact_for_pat(A, A_pat)

% X = ssa_compute_exact_for_pat(A, A_pat)
%
% Wrapper function over all the sparse spectral approximation (ssa)
% functionality.  This one uses the "exact" algorithm which means
% no binning is done, and null-space constraints are imposed while
% minimization rather than in the last step.
%
% Given
%
%   a matrix A
%   a pattern A_pat of same size as A with zeros and ones
%
% this function computes a matrix X that is sparse (like A_pat is)
% and spectrally close to A, especially in the lower end of the singular
% value spectrum.  Even if A is complex, A_pat must be real.
%
% No error checking is performed in this function.  It is assumed that
% the functions that are called will do the checks.

[pinv_A, rnull, lnull] = pinv_rrqr(A);

pinv_ATA = pinv_A * pinv_A';
pinv_AAT = pinv_A' * pinv_A;

X = ssa_minimization_with_action(A, sparse(A_pat), pinv_ATA, pinv_AAT, rnull, lnull);
