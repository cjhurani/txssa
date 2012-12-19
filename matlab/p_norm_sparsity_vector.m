function v_pat = p_norm_sparsity_vector(v, ratio, p, varargin)

% v_pat = p_norm_sparsity_vector(v, ratio, p [, min_num_nnz])
%
% Computes a pattern vector v_pat which is a p-norm sparsity pattern for v.
%
% Whether v is real or complex, v_pat is sparse double real of the same shape
% and contains only 0 or 1.
%
% Pre-conditions: ratio should be in [0,1].  p should be in [0,inf]
%
% varargin can be min_num_nnz.
% min_num_nnz is the minimum number of non-zeros needed in the pattern.
% Default is 0.

optargin = size(varargin, 2);
assert(optargin <= 1, 'p_norm_sparsity_vector: Too many optional input arguments (%d)', optargin);

assert(isvector(v), 'p_norm_sparsity_vector: first input (v) must be a 1-D vector');

assert(numel(ratio) == 1,        'p_norm_sparsity_vector: second input (ratio) must be a number');
assert(isreal(ratio),            'p_norm_sparsity_vector: second input (ratio) must be real');
assert(isfloat(ratio),           'p_norm_sparsity_vector: second input (ratio) must be a floating point number', ratio);
assert(0 <= ratio && ratio <= 1, 'p_norm_sparsity_vector: second input (ratio = %f) must be in [0,1]', ratio);

assert(numel(p) == 1, 'p_norm_sparsity_vector: third input must be a number');
assert(isreal(p),     'p_norm_sparsity_vector: third input must be real');
assert(0 <= p,        'p_norm_sparsity_vector: third input (p = %f) must not be less than 0', p);

v = abs(v);
nnz_v = sum(0 < v);

if(nnz_v == 0)
    v_pat = sparse(size(v,1),size(v,2));
    return;
end

if(optargin == 0 || isempty(varargin{1}))
    min_num_nnz = 0;
else
    min_num_nnz = varargin{1};

    if(numel(min_num_nnz) ~= 1 || ~isreal(min_num_nnz) || int32(min_num_nnz) ~= min_num_nnz || nnz_v < min_num_nnz || min_num_nnz < 0)
        assert(false, 'p_norm_sparsity_vector: min_num_nnz should numerically be an integer within range')
    end
end

if(0 < p && p < Inf)
    v = v .^ p;
end

% Ascending sort to match C++ version.  Otherwise if multiple entries
% are equal, the one that actually is used depends on column id or row
% id.
dim = 2;
if(size(v,2) < size(v,1))
    dim = 1;
end
[v_vals j] = sort(v, dim, 'ascend');

n = length(v);

if(p == 0)
    c = ceil(ratio * nnz_v);
    % if ratio * nnz_v is an integer + "something positive very small",
    % then reduce c because ceil would have increased it.
    if(c - ratio * nnz_v > 1 - 1000*eps)  % MAGIC CONSTANT
        c = c - 1;
    end
    num_needed = max(min_num_nnz, c);
elseif(p == Inf)
    if(min_num_nnz == nnz_v)
        num_needed = min_num_nnz;
    else
        v_norm = v_vals(end);
        discarded_bound = (1-ratio)*v_norm;
        num_needed = n - upper_bound(v_vals(1:end-min_num_nnz), discarded_bound) + 1;
    end
else
    csum = cumsum(v_vals);
    sum_threshold = csum(end) * (1 - ratio) ^ max(1,p);
    ub = upper_bound(csum, sum_threshold); % a(1:ub-1) <= sum_threshold 
    % ub will be in [1,n]
    num_needed = max(min_num_nnz, n - ub + 1);
end

v_pat = sparse([],[],[],size(v,1),size(v,2));
v_pat(j(n - num_needed + 1:n)) = 1;
