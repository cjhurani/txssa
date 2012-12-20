function A_id = bin_sparse_matrix(A, A_pat, max_num_bins)

% A_id = bin_sparse_matrix(A, A_pat, max_num_bins)
% 
% Computes A_id which is a sparse matrix and of class double that contains
% the one-based DOF ids corresponding to A and A_pat.  The DOF ids can be binned
% depending on the value of max_num_bins.  Binned means more than one non-zero
% entry can have the same value.
%
% 'A_pat' must be a real square/rectangular sparse matrix with zeros and
% ones.  A can be single/double real/complex matrix of the same size
% as 'A_pat'.
%
% A is the matrix whose entries are used for binning.  If max_num_bins is 0,
% then values in A are not referenced and can be anything.  Only 'A_pat' is
% used then.  We do use the field of A (real or complex).
%
% max_num_bins must be non-negative.  If max_num_bins == 0, no binning is
% performed.

assert(issparse(A_pat), 'bin_sparse_matrix: second input (A_pat) must be sparse');
assert(isreal(A_pat),   'bin_sparse_matrix: second input (A_pat) must be real');

assert(~isempty(max_num_bins), 'bin_sparse_matrix: third input (max_num_bins) must have some elements');
assert(isreal(max_num_bins),   'bin_sparse_matrix: third input (max_num_bins) must be real');
assert(0 <= max_num_bins,      'bin_sparse_matrix: third input (max_num_bins) must be non-negative');
assert(max_num_bins == int32(max_num_bins), 'bin_sparse_matrix: third input (max_num_bins) must be integral');

[n1 n2] = size(A_pat);

if(0 < max_num_bins && (size(A,1) ~= n1 || size(A,2) ~= n2))
    assert(false, 'bin_sparse_matrix: If 0 < max_num_bins, A must be the same size as A_pat')
end

[A_i A_j] = find(A_pat);
A_pat_nnz = nnz(A_pat);

is_real_A = isreal(A);

perturb_fuzz = 1E3;
separated_at = 0;

if(is_real_A)
    if(max_num_bins == 0)
        A_id = sparse(A_i, A_j, 1:A_pat_nnz, n1, n2);
    else
        nz_A = zeros(A_pat_nnz,1);

        for k = 1:A_pat_nnz
            nz_A(k) = A(A_i(k), A_j(k));
        end

        [min_left, max_left, min_right, max_right] = separated_min_max(nz_A, separated_at, perturb_fuzz);
        bin_ids = binned_with_separated_min_max(nz_A, max_num_bins, min_left, max_left, min_right, max_right, separated_at);
        new_ids = bin_mapping(bin_ids);
        A_id = sparse(A_i, A_j, double(new_ids), n1, n2);
    end
else
    if(max_num_bins == 0)
        real_A_id = sparse(A_i, A_j, 1:A_pat_nnz, n1, n2);
        imag_A_id = sparse(A_i, A_j, (A_pat_nnz+1):(2*A_pat_nnz), n1, n2);
    else
        real_A = real(A);
        imag_A = imag(A);

        real_nz_A = zeros(A_pat_nnz,1);

        for k = 1:A_pat_nnz
            real_nz_A(k) = real_A(A_i(k), A_j(k));
        end

        imag_nz_A = zeros(A_pat_nnz,1);

        for k = 1:A_pat_nnz
            imag_nz_A(k) = imag_A(A_i(k), A_j(k));
        end

        % bin real and imaginary separately

        [min_left, max_left, min_right, max_right] = separated_min_max(real_nz_A, separated_at, perturb_fuzz);
        bin_ids = binned_with_separated_min_max(real_nz_A, max_num_bins, min_left, max_left, min_right, max_right, separated_at);
        new_ids = bin_mapping(bin_ids);
        max_real_id = max(new_ids);
        real_A_id = sparse(A_i, A_j, double(new_ids), n1, n2);

        [min_left, max_left, min_right, max_right] = separated_min_max(imag_nz_A, separated_at, perturb_fuzz);
        bin_ids = binned_with_separated_min_max(imag_nz_A, max_num_bins, min_left, max_left, min_right, max_right, separated_at);
        new_ids = bin_mapping(bin_ids);
        imag_A_id = sparse(A_i, A_j, double(new_ids + max_real_id), n1, n2);
    end

    % complex(a,b) doesn't allow sparse a and b, so use i.
    A_id = real_A_id + 1i * imag_A_id;
end
