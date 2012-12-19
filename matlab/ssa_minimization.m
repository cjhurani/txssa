function X = ssa_minimization(A, A_id, varargin)

% X = ssa_minimization(A, A_id [, pinv_ATA] [, pinv_AAT])
%
% 'ssa' stands for sparse spectral approximation.
%
% Computes a sparse X of the same size and field as A such that X approximates
% the lower spectrum of A well.  All arguments must be either real or all must be
% complex.  Null-spaces of A are not imposed on X here.  This is because of
% binning (whether it is actually done or not).
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
% If not provided, pinv_ATA and pinv_AAT are computed inside the function.
%
% If pinv_ATA and pinv_AAT are computed from A or are given based on A, the
% optimization process minimizes
%
%   (1/2) (|| (X - A) pA ||_F^2 + || pA (X - A) ||_F^2)
%
% where pA = pinv(A), subject to sparsity and binning constraints on X specified
% in A_id.

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

optargin = size(varargin, 2);

assert(optargin <= 2, 'ssa_minimization: Too many optional input arguments (%d)', optargin);

compute_pinv_AAT = false;
compute_pinv_ATA = false;

[n1 n2] = size(A);

if(optargin <= 1 || isempty(varargin{2}))
    compute_pinv_AAT = true;
else
    pinv_AAT = varargin{2};
    assert(all(size(pinv_AAT) == n1) && ~xor(is_real_A, isreal(pinv_AAT)), 'ssa_minimization: pinv_AAT is not of correct size or field');
end

if(optargin == 0 || isempty(varargin{1}))
    compute_pinv_ATA = true;
else
    pinv_ATA = varargin{1};
    assert(all(size(pinv_ATA) == n2) && ~xor(is_real_A, isreal(pinv_ATA)), 'ssa_minimization: pinv_ATA is not of correct size or field');
end

if(compute_pinv_AAT || compute_pinv_ATA)
    pinv_A = pinv_rrqr(A);

    if(compute_pinv_AAT)
        pinv_AAT = pinv_A' * pinv_A;
    end

    if(compute_pinv_ATA)
        pinv_ATA = pinv_A * pinv_A';
    end
end

[LS_A, LS_b] = ssa_system_no_null(A, A_id, pinv_ATA, pinv_AAT);

if(size(LS_A,1) ~= 0)
    LS_x = pinv_rrqr(LS_A) * LS_b;
else
    LS_x = [];
end

X = ssa_unknown_to_matrix(LS_x, A_id);

