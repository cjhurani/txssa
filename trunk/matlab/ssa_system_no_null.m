function [LS_A, LS_b] = ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT)

% [LS_A, LS_b] = ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT)
%
% 'ssa' stands for sparse spectral approximation.
%
% Computes the system matrix and system rhs for the minimization problem
% without any null-space constraints.  LS_A and LS_b are real even if
% A is complex.  All arguments must be either real or all must be complex.
% 
% A_id is a sparse, double matrix of size(A) to provide sparsity and binning
% information.  It must have non-negative integer entries where 0 means location
% of 0, positive values mean id of unknowns.  If complex, these conditions apply
% to the real and imaginary parts separately.  The two parts can have
% overlapping ids too.
%
% pinv_ATA = pinv(A'*A) and pinv_AAT = pinv(A*A').  We don't check that
% these arguments are actually related to A.  So in practice one can
% pass something else.  One reason we don't compute them inside using pinv
% is because the caller may compute them more efficiently outside if
% more is known about A or because pinv of an SPSD matrix is needed.
%
% If pinv_ATA and pinv_AAT are computed from A or are given based on A, the
% system matrix and system rhs correspond to the first-order optimality
% conditions for the following misfit
%
%   (1/2) (|| (X - A) pA ||_F^2 + || pA (X - A) ||_F^2)
%
% where pA = pinv(A), subject to sparsity and binning constraints on X specified
% in A_id.

is_real_A = isreal(A);

assert(nnz(A_id) == 0 || ~xor(is_real_A, isreal(A_id)), 'ssa_system_no_null: A_id must be a of the same field as A');
assert(all(size(A_id) == size(A)), 'ssa_system_no_null: A_id must be of the size as A');
assert(issparse(A_id), 'ssa_system_no_null: A_id must be sparse');

if(is_real_A)
    assert(0 <= min(min(A_id)), 'ssa_system_no_null: A_id contains negative integers');
else
    real_A_id = real(A_id);
    imag_A_id = imag(A_id);

    assert(0 <= min(min(real_A_id)) && 0 <= min(min(imag_A_id)), 'ssa_system_no_null: A_id contains negative integers');
end

[n1 n2] = size(A);

assert(all(size(pinv_AAT) == n1) && ~xor(is_real_A, isreal(pinv_AAT)), 'ssa_system_no_null: pinv_AAT is not of correct size or field');
assert(all(size(pinv_ATA) == n2) && ~xor(is_real_A, isreal(pinv_ATA)), 'ssa_system_no_null: pinv_ATA is not of correct size or field');

if(is_real_A)
    [A_id_i A_id_j A_id_v] = find(A_id);
else
    real_pinv_AAT = real(pinv_AAT);
    imag_pinv_AAT = imag(pinv_AAT);
    real_pinv_ATA = real(pinv_ATA);
    imag_pinv_ATA = imag(pinv_ATA);

    [real_A_id_i real_A_id_j real_A_id_v] = find(real_A_id);
    [imag_A_id_i imag_A_id_j imag_A_id_v] = find(imag_A_id);
    A_id_v = [real_A_id_v; imag_A_id_v];
end

A_id_uniq = unique(A_id_v);
N = length(A_id_uniq);

% Check that dof ids begin from 1 and have no 'holes'.
% We assume that all entries are integers.
% unique() would have returned a sorted array.
if(~isempty(A_id_uniq))
    assert(~(N == 0 || A_id_uniq(1) ~= 1 || A_id_uniq(end) ~= N), 'ssa_system_no_null: Non-zero entries of A_id not compatible');
end

% number of most frequently occuring dof id
max_id_occurences = sum(A_id_v == mode(A_id_v));

% A_id_loc_row(id, occurence) contains the row of a location (i,j)
% A_id_loc_col(id, occurence) contains the col of a location (i,j)
% A_id_loc_prt(id, occurence) contains the real/imag location. 1 = real, 2 = imag.
% A_id_count(id) contains the number of occurences.

A_id_loc_row = zeros(N, max_id_occurences);
A_id_loc_col = zeros(N, max_id_occurences);
if(~is_real_A)
    A_id_loc_prt = zeros(N, max_id_occurences, 'uint8');
end
A_id_loc_next = ones(N, 1);

if(is_real_A)
    A_id_nnz = nnz(A_id);

    for k = 1:A_id_nnz
        A_id_val = A_id_v(k);
        curr = A_id_loc_next(A_id_val);
        A_id_loc_row(A_id_val, curr) = A_id_i(k);
        A_id_loc_col(A_id_val, curr) = A_id_j(k);
        A_id_loc_next(A_id_val) = curr + 1;
    end
else
    real_A_id_nnz = nnz(real_A_id);

    for k = 1:real_A_id_nnz
        real_A_id_val = real_A_id_v(k);
        curr = A_id_loc_next(real_A_id_val);
        A_id_loc_row(real_A_id_val, curr) = real_A_id_i(k);
        A_id_loc_col(real_A_id_val, curr) = real_A_id_j(k);
        A_id_loc_prt(real_A_id_val, curr) = 1;
        A_id_loc_next(real_A_id_val) = curr + 1;
    end

    imag_A_id_nnz = nnz(imag_A_id);

    for k = 1:imag_A_id_nnz
        imag_A_id_val = imag_A_id_v(k);
        curr = A_id_loc_next(imag_A_id_val);
        A_id_loc_row(imag_A_id_val, curr) = imag_A_id_i(k);
        A_id_loc_col(imag_A_id_val, curr) = imag_A_id_j(k);
        A_id_loc_prt(imag_A_id_val, curr) = 2;
        A_id_loc_next(imag_A_id_val) = curr + 1;
    end
end

A_id_count = A_id_loc_next - 1;

A__pinv_ATA = A * pinv_ATA;
pinv_AAT__A = pinv_AAT * A;

if(~is_real_A)
    real_A__pinv_ATA = real(A__pinv_ATA);
    imag_A__pinv_ATA = imag(A__pinv_ATA);
    real_pinv_AAT__A = real(pinv_AAT__A);
    imag_pinv_AAT__A = imag(pinv_AAT__A);
end

% minimization system matrix and rhs vector

class_A = class(A);

LS_A = zeros(N, class_A);
LS_b = zeros(N, 1, class_A);

% form LS system

% d2J11 = kron(real_pinv_ATA, eye(n1)) + kron(eye(n2), real_pinv_AAT)
% d2J12 = kron(imag_pinv_ATA, eye(n1)) - kron(eye(n2), imag_pinv_AAT)
% for unconstrained complex problem d2J = [d2J11 d2J12; -d2J12 d2J11];
% for unconstrained real    problem d2J = d2J11;

for i = 1:N
    for j = 1:N

        tmp = 0;
        for K = 1:A_id_count(i)
            for L = 1:A_id_count(j)

                u_row = A_id_loc_row(i, K);
                u_col = A_id_loc_col(i, K);

                v_row = A_id_loc_row(j, L);
                v_col = A_id_loc_col(j, L);

                if(is_real_A)
                    if(u_col == v_col)
                        tmp = tmp + pinv_AAT(u_row, v_row);
                    end

                    if(u_row == v_row)
                        tmp = tmp + pinv_ATA(u_col, v_col);
                    end
                else
                    u_prt = A_id_loc_prt(i, K);
                    v_prt = A_id_loc_prt(j, L);
                    
                    if(u_col == v_col)
                        if(    u_prt == 1 && v_prt == 1)
                            tmp = tmp + real_pinv_AAT(u_row, v_row);
                        elseif(u_prt == 1 && v_prt == 2)
                            tmp = tmp - imag_pinv_AAT(u_row, v_row);
                        elseif(u_prt == 2 && v_prt == 1)
                            tmp = tmp + imag_pinv_AAT(u_row, v_row);
                        elseif(u_prt == 2 && v_prt == 2)
                            tmp = tmp + real_pinv_AAT(u_row, v_row);
                        end
                    end

                    if(u_row == v_row)
                        if(    u_prt == 1 && v_prt == 1)
                            tmp = tmp + real_pinv_ATA(u_col, v_col);
                        elseif(u_prt == 1 && v_prt == 2)
                            tmp = tmp + imag_pinv_ATA(u_col, v_col);
                        elseif(u_prt == 2 && v_prt == 1)
                            tmp = tmp - imag_pinv_ATA(u_col, v_col);
                        elseif(u_prt == 2 && v_prt == 2)
                            tmp = tmp + real_pinv_ATA(u_col, v_col);
                        end
                    end
                end
            end
        end

        LS_A(i,j) = tmp;
    end
end

%for fully unconstrained complex problem
% rhs = [
% real_A__pinv_ATA(:) + real_pinv_AAT__A(:)
% imag_A__pinv_ATA(:) + imag_pinv_AAT__A(:)
% ]

%for fully unconstrained real problem
% rhs = [
% real_A__pinv_ATA(:) + real_pinv_AAT__A(:)
% ]

for i = 1:N
    tmp = 0;
    for K = 1:A_id_count(i)
        row = A_id_loc_row(i, K);
        col = A_id_loc_col(i, K);

        if(is_real_A)
            tmp = tmp + A__pinv_ATA(row, col) + pinv_AAT__A(row, col);
        else
            prt = A_id_loc_prt(i, K);
            if(prt == 1)
                tmp = tmp + real_A__pinv_ATA(row, col) + real_pinv_AAT__A(row, col);
            elseif(prt == 2)
                tmp = tmp + imag_A__pinv_ATA(row, col) + imag_pinv_AAT__A(row, col);
            else
                assert(false, 'ssa_system_no_null: Internal error - Bad value in A_id_loc_prt');
            end
        end
    end
    LS_b(i) = tmp;
end

