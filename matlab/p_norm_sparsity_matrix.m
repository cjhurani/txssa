function A_pat = p_norm_sparsity_matrix(A, ratio, p, varargin)

% A_pat = p_norm_sparsity_matrix(A, ratio, p [, min_per_row] [, min_per_col])
%
% Computes a pattern matrix A_pat which is a p-norm sparsity pattern for A.
%
% Whether A is real or complex, A_pat is sparse double real or same shape and
% contains only 0 or 1.
%
% ratio should be in [0,1].  p should be in [0,inf]
%
% varargin can be min_per_row and min_per_col.
% min_per_row is the minimum number of non-zeros needed per row.
% min_per_col is the minimum number of non-zeros needed per column.
% Defaults are 0.

optargin = size(varargin, 2);

assert(optargin <= 2, 'p_norm_sparsity_matrix: Too many optional input arguments (%d)', optargin);

assert(numel(ratio) == 1, 'p_norm_sparsity_matrix: second input (ratio) must be a number');
assert(isreal(ratio),     'p_norm_sparsity_matrix: second input (ratio) must be real');
assert(isfloat(ratio),    'p_norm_sparsity_matrix: second input (ratio) must be a floating point number', ratio);
assert(0 <= ratio && ratio <= 1, 'p_norm_sparsity_matrix: second input (ratio = %f) must be in [0,1]', ratio);

assert(numel(p) == 1, 'p_norm_sparsity_matrix: third input must be a number');
assert(isreal(p),     'p_norm_sparsity_matrix: third input must be real');
assert(0 <= p,        'p_norm_sparsity_matrix: third input (p = %f) must not be less than 0', p);

if(optargin == 0 || isempty(varargin{1}))
    min_per_row = 0;
else
    min_per_row = varargin{1};
end

if(optargin <= 1 || isempty(varargin{2}))
    min_per_col = 0;
else
    min_per_col = varargin{2};
end

[n1 n2] = size(A);

if(numel(min_per_row) ~= 1 || ~isreal(min_per_row) || ...
    int32(min_per_row) ~= min_per_row || n2 < min_per_row || min_per_row < 0)
    assert(false, 'p_norm_sparsity_matrix: min_per_row should numerically be an integer within range')
end

if(numel(min_per_col) ~= 1 || ~isreal(min_per_col) || ...
    int32(min_per_col) ~= min_per_col || n1 < min_per_col || min_per_col < 0)
    assert(false, 'p_norm_sparsity_matrix: min_per_col should numerically be an integer within range')
end

% For square matrices, check if the absolute value is almost symmetric, and if
% so, make the matrix exactly symmetric.  This is so that the sparsity pattern
% is symmetric too for such matrices.  Note that the code below automatically
% takes care of real-symmetric, complex-symmetric, complex-Hermitian, and their
% skew counterparts.  That's because in each case, the entry-wise absolute value
% is symmetric.

if(n1 == n2)
    A = abs(A);

    sym_diff = max(max(abs(A - A')));
    a_max    = max(abs(A(:)));

    class_a = class(A);
    make_sym = false;

    if(strcmp(class_a, 'single') && sym_diff <= 100 * eps(class_a) * a_max)  % MAGIC CONSTANT
        make_sym = true;
    elseif(strcmp(class_a, 'double') && sym_diff <= 1000 * eps(class_a) * a_max)  % MAGIC CONSTANT
        make_sym = true;
    end

    if(make_sym)
        A = 0.5*(A + A');

        min_per_row = max(min_per_row, min_per_col);
        min_per_col = min_per_row;
    end
end

A_pat = sparse([],[],[],n1,n2);

for i = 1:n1
    i_row_pat = p_norm_sparsity_vector(A(i,:), ratio, p, min_per_row);
    A_pat(i,:) = i_row_pat;
end

for j = 1:n2
    j_col_pat = p_norm_sparsity_vector(A(:,j), ratio, p, min_per_col);
    A_pat(:,j) = max(A_pat(:,j), j_col_pat);
end
