function X = ssa_minimization_with_action(A, A_pat, pinv_ATA, pinv_AAT, R_mat, L_mat)

% X = ssa_minimization_with_action(A, A_pat, pinv_ATA, pinv_AAT, R_mat, L_mat)
%
% 'ssa' stands for sparse spectral approximation.
%
% Computes a sparse X of the same size and field as A such that X approximates
% the lower spectrum of A well.  Whether A is real or complex, A_pat must be
% real with zeros and ones.
%
% pinv_ATA = pinv(A'*A) and pinv_AAT = pinv(A*A').  One reason we don't compute
% them inside using pinv is because the caller may compute them more efficiently
% outside if more is known about A or because pinv of an SPSD matrix is needed.
% Another reason is that null-spaces of A are used in computing the pattern
% (A_pat) and they are by-product of computing the pseudoinverse.
%
% The computed X is such that (X - A) * R_mat = 0, and (X' - A') * L_mat = 0.
% Thus, R_mat and L_mat must be of appropriate sizes.  They could be
% empty too.  Typically, R_mat and L_mat form a basis for the right and left
% null-spaces of A, respectively.  The optimization process minimizes
%
% The optimization process minimizes
%
%   (1/2) (|| (X - A) pA ||_F^2 + || pA (X - A) ||_F^2)
%
% where pA = pinv(A), subject to sparsity constraints on X specified in A_pat
% and matrix-action constraints coming from A, R_mat, and L_mat.

assert(isreal(A_pat), 'ssa_minimization_with_action: A_pat must be real');
assert(issparse(A_pat), 'ssa_minimization_with_action: A_pat must be sparse');
assert(all(size(A_pat) == size(A)), 'ssa_minimization_with_action: A_pat must be of the size as A');

[n1 n2] = size(A);

unique_A_pat = unique(A_pat(:));
unique_A_pat(unique_A_pat == 0) = [];

assert(isempty(unique_A_pat) || (numel(unique_A_pat) == 1 && unique_A_pat == 1), 'ssa_minimization_with_action: A_pat is not made of zeros and ones only');

assert(all(size(pinv_AAT) == n1) && ~xor(isreal(A), isreal(pinv_AAT)), 'ssa_minimization_with_action: pinv_AAT is not of correct size or field');
assert(all(size(pinv_ATA) == n2) && ~xor(isreal(A), isreal(pinv_ATA)), 'ssa_minimization_with_action: pinv_ATA is not of correct size or field');

[nR1 nR2] = size(R_mat);
[nL1 nL2] = size(L_mat);

assert(n2 == nR1, 'ssa_minimization_with_action: A * R_mat not possible');
assert(n1 == nL1, 'ssa_minimization_with_action: ctranspose(A) * L_mat not possible');
assert(nR2 == 0 || ~xor(isreal(A), isreal(R_mat)), 'ssa_minimization_with_action: R_mat must be a of the same field as A');
assert(nL2 == 0 || ~xor(isreal(A), isreal(L_mat)), 'ssa_minimization_with_action: L_mat must be a of the same field as A');

max_num_bins = 0;
A_id = bin_sparse_matrix(A, A_pat, max_num_bins);

[LS_A, LS_b] = ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT);
[C_lhs, C_rhs] = ssa_action_constraint(A, A_id, R_mat, L_mat);

if(~isequal(class(LS_A), 'single'))
    LS_A = sparse(LS_A);
    LS_full_lhs = [LS_A C_lhs'; C_lhs sparse(length(C_rhs), length(C_rhs))];
else
    LS_full_lhs = [LS_A C_lhs'; C_lhs zeros(length(C_rhs), length(C_rhs))];
end

LS_full_rhs = [LS_b; C_rhs];

if(size(LS_full_lhs,1) ~= 0)
    LS_x = pinv_rrqr(LS_full_lhs) * LS_full_rhs;
else
    LS_x = [];
end

X = ssa_unknown_to_matrix(LS_x(1:length(LS_b)), A_id);
