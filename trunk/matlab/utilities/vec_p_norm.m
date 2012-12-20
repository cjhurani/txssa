function val = vec_p_norm(v, p)

% For p in [0,inf] and real or complex vector v, vec_p_norm(v, p)
% computes the p-"norm" of v.  The definition of p-"norm":
%
% ||v||_p = p-"norm" of a collection of numbers [v_1, v_2, ..., v_n] =
% (sum_{i=1}^n |v_i|^p)^(1/p)  for 1 <= p < inf
%  sum_{i=1}^n |v_i|^p         for 0 <  p < 1    Note: no 1/p.
%  max_{i} |v_i|               for p = inf
%  number of non-zero |v_i|    for p = 0
%
% It is a "norm" only when 1 <= p <= inf.  Hence, for a general p, it
% should be called a p-"norm".

assert(min(size(v)) <= 1, 'vec_p_norm: first input must be 0 or 1-D');
assert(numel(p) == 1,     'vec_p_norm: second input must be a single element');
assert(isreal(p),         'vec_p_norm: second input must be real');
assert(0 <= p,            'vec_p_norm: second input (%f) must not be less than 0', p);

if(p == 0)
    val = nnz(v);
elseif(p == 1)
    val = sum(abs(v));
elseif(p == Inf)
    val = max(abs(v));
else
    val = sum(abs(v).^p)^min(1, 1/p);
end
