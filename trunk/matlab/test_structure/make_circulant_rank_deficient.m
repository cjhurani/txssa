function X = make_circulant_rank_deficient(A, even_rank_to_remove)

% X = make_circulant_rank_deficient(A, even_rank_to_remove)
%
% Output a circulant or skew_circulant real matrix that is near
% the input matrix and has "even_rank_to_remove" less rank.
% even_rank_to_remove must be even.

if(ndims(A) ~= 2 || ~isreal(A))
    assert(false, 'make_circulant_rank_deficient: Input A is not 2-D or not real');
end

assert(numel(even_rank_to_remove) == 1,     'make_circulant_rank_deficient: even_rank_to_remove must be a single element');
assert(isreal(even_rank_to_remove),         'make_circulant_rank_deficient: even_rank_to_remove must be real');
assert(0 <= even_rank_to_remove,            'make_circulant_rank_deficient: even_rank_to_remove not be negative');
assert(even_rank_to_remove <= min(size(A)), 'make_circulant_rank_deficient: even_rank_to_remove is larger than minimum size of A');
assert(mod(even_rank_to_remove,2) == 0,     'make_circulant_rank_deficient: even_rank_to_remove is not even');

if(min(size(A)) == 0 && even_rank_to_remove == 0)
    X = full(A);
    return;
end

[u s v] = svd(full(A));

diag_s = diag(s);

% for real even sized circulant, an example pattern for singular values is
% count     value
%     1     0.4799
%     1     1.9091
%     2     0.2580
%     2     0.8433
%...
%     2     1.0175

% for real odd sized circulant and real odd skew circulant, an example
% pattern for singular values is
% count     value
%     1     0.1976
%     2     0.3193
%     2     0.6418
%     2     0.8433
%...
%     2     1.0175


if(diag_s(1) == 0)
    if(even_rank_to_remove == 0)
        X = full(A);
        return;
    else
        assert(false, 'make_circulant_rank_deficient: even_rank_to_remove is larger than rank(A)');
    end
end

tol = min(size(A)) * diag_s(1) * eps(class(A));
rnk = sum(tol < diag_s);

assert(even_rank_to_remove <= rnk, 'make_circulant_rank_deficient: even_rank_to_remove is larger than rank(A) = %d', rnk);

diff_diag_s = diff(diag_s);

diff_diag_s(abs(diff_diag_s) < tol) = 0;

for i = 1:length(diff_diag_s)
    if(diff_diag_s(i) == 0 && even_rank_to_remove > 0)
        diag_s(i) = 0;
        diag_s(i+1) = 0;
        even_rank_to_remove = even_rank_to_remove - 2;
    end
end

if(even_rank_to_remove > 0)
    assert(false, 'make_circulant_rank_deficient: Could not introduce rank deficiency');
end

s = diag(diag_s);

X = u*s*v';
