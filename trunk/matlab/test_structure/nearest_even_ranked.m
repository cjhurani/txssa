function X = nearest_even_ranked(A, even_rank)

% X = nearest_even_ranked(A, even_rank)
%
% Output the nearest matrix to A with rank 'even_rank'.
% even_rank must be even and <= rank(A).

if(ndims(A) ~= 2)
    assert(false, 'nearest_even_ranked: Input A is not 2-D');
end

assert(numel(even_rank) == 1,     'nearest_even_ranked: even_rank must be a single element');
assert(isreal(even_rank),         'nearest_even_ranked: even_rank must be real');
assert(0 <= even_rank,            'nearest_even_ranked: even_rank not be negative');
assert(even_rank <= min(size(A)), 'nearest_even_ranked: even_rank is larger than minimum size of A');
assert(mod(even_rank,2) == 0,     'nearest_even_ranked: even_rank is not even');

if(min(size(A)) == 0 && even_rank == 0)
    X = full(A);
    return;
end

[u s v] = svd(full(A));

diag_s = diag(s);

if(diag_s(1) == 0)
    if(even_rank == 0)
        X = full(A);
        return;
    else
        assert(false, 'nearest_even_ranked: even_rank is larger than rank(A)');
    end
end

tol = min(size(A)) * diag_s(1) * eps(class(A));
rnk = sum(tol < diag_s);

assert(even_rank <= rnk, 'nearest_even_ranked: even_rank is larger than rank(A) = %d', rnk);

s(even_rank+1:end, even_rank+1:end) = 0;

X = u*s*v';
