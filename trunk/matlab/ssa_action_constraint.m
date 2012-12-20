function [C_lhs, C_rhs] = ssa_action_constraint(A, A_id, R_mat, L_mat)

% [C_lhs, C_rhs] = ssa_action_constraint(A, A_id, R_mat, L_mat)
%
% Computes a sparse C_lhs and a vector C_rhs of the same field as A such
% that C_lhs and C_rhs represent the constraints to be imposed for
% (X - A)*R_mat = 0 and (X - A)'*L_mat = 0 for some unknown matrix
% X on which such constraints are to be imposed.
%
% A_id is a sparse, double matrix of size(A) to provide sparsity and binning
% information.  It must have non-negative integer entries where 0 means location
% of 0, positive values mean id of unknowns.  If complex, these conditions apply
% to the real and imaginary parts separately.  The two parts can have
% overlapping ids too.

is_real_A = isreal(A);
assert(nnz(A_id) == 0 || ~xor(is_real_A, isreal(A_id)), 'ssa_minimization: A_id must be a of the same field as A');
assert(all(size(A_id) == size(A)), 'ssa_minimization: A_id must be of the size as A');
assert(issparse(A_id), 'ssa_minimization: A_id must be sparse');

if(is_real_A)
    assert(0 <= min(min(A_id)), 'ssa_minimization: A_id contains negative integers');
else
    real_A_id = real(A_id);
    imag_A_id = imag(A_id);

    assert(0 <= min(min(real_A_id)) && 0 <= min(min(imag_A_id)), 'ssa_minimization: A_id contains negative integers');
end

[nR1 nR2] = size(R_mat);
[nL1 nL2] = size(L_mat);
[n1 n2] = size(A);

assert(n2 == nR1, 'ssa_action_constraint: A * R_mat not possible');
assert(n1 == nL1, 'ssa_action_constraint: ctranspose(A) * L_mat not possible');
assert(nR2 == 0 || ~xor(isreal(A), isreal(R_mat)), 'ssa_action_constraint: R_mat must be a of the same field as A');
assert(nL2 == 0 || ~xor(isreal(A), isreal(L_mat)), 'ssa_action_constraint: L_mat must be a of the same field as A');

if(is_real_A)
    [tmp tmp A_id_v] = find(A_id);
else
    [tmp tmp real_A_id_v] = find(real_A_id);
    [tmp tmp imag_A_id_v] = find(imag_A_id);
    A_id_v = [real_A_id_v; imag_A_id_v];
end

A_id_uniq = unique(A_id_v);
N = length(A_id_uniq);

A_R_mat = A * R_mat;
AT_L_mat = A' * L_mat;

if(is_real_A)
    C_lhs = sparse(n1 * nR2 + n2 * nL2, N);
    C_rhs = zeros(size(C_lhs,1),1);

    for i = 1:n1
        [tmp col_ids dof_ids] = find(A_id(i,:));
        for jj = 1:length(col_ids)
            j = col_ids(jj);
            dof = dof_ids(jj);
            for k = 1:nR2
                C_row = i + (k-1)*n1;
                C_lhs(C_row, dof) = R_mat(j,k);
                C_rhs(C_row) = A_R_mat(i,k);
            end
        end
    end

    for j = 1:n2
        [row_ids tmp dof_ids] = find(A_id(:,j));
        for ii = 1:length(row_ids)
            i = row_ids(ii);
            dof = dof_ids(ii);
            for k = 1:nL2
                C_row = n1 * nR2 + j + (k-1)*n2;
                C_lhs(C_row, dof) = L_mat(i,k);
                C_rhs(C_row) = AT_L_mat(j,k);
            end
        end
    end
else

    C_lhs = sparse(2*(n1 * nR2 + n2 * nL2), N);
    C_rhs = zeros(size(C_lhs,1),1);

    for i = 1:n1
        [tmp col_ids dof_ids] = find(A_id(i,:));
        for jj = 1:length(col_ids)
            j = col_ids(jj);
            real_dof = real(dof_ids(jj));
            imag_dof = imag(dof_ids(jj));
            for k = 1:nR2
                C_row_real = 2*(i + (k-1)*n1) - 1;
                C_row_imag = C_row_real + 1;
                C_lhs(C_row_real, real_dof) = C_lhs(C_row_real, real_dof) + real(R_mat(j,k));
                C_lhs(C_row_real, imag_dof) = C_lhs(C_row_real, imag_dof) - imag(R_mat(j,k));
                C_lhs(C_row_imag, real_dof) = C_lhs(C_row_imag, real_dof) + imag(R_mat(j,k));
                C_lhs(C_row_imag, imag_dof) = C_lhs(C_row_imag, imag_dof) + real(R_mat(j,k));
                C_rhs(C_row_real) = real(A_R_mat(i,k));
                C_rhs(C_row_imag) = imag(A_R_mat(i,k));
            end
        end
    end

    for j = 1:n2
        [row_ids tmp dof_ids] = find(A_id(:,j));
        for ii = 1:length(row_ids)
            i = row_ids(ii);
            real_dof = real(dof_ids(ii));
            imag_dof = imag(dof_ids(ii));
            for k = 1:nL2
                C_row_real = 2*(n1 * nR2 + j + (k-1)*n1) - 1;
                C_row_imag = C_row_real + 1;
                C_lhs(C_row_real, real_dof) = C_lhs(C_row_real, real_dof) + real(L_mat(i,k));
                C_lhs(C_row_real, imag_dof) = C_lhs(C_row_real, imag_dof) + imag(L_mat(i,k));
                C_lhs(C_row_imag, real_dof) = C_lhs(C_row_imag, real_dof) - imag(L_mat(i,k));
                C_lhs(C_row_imag, imag_dof) = C_lhs(C_row_imag, imag_dof) + real(L_mat(i,k));
                C_rhs(C_row_real) = real(AT_L_mat(j,k));
                C_rhs(C_row_imag) = imag(AT_L_mat(j,k));
            end
        end
    end
end
