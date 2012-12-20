function diff = diff_skew_circulant(A)

% diff = diff_skew_circulant(A)
%
% Checks how far is a matrix A different from being skew_circulant.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2))
    assert(false, 'diff_circulant: Input A is not 2-D or not square');
end

if(any(size(A) == 0))
    diff = 0;
    return;
end

n = size(A,1);

z = sparse(n-1,1);
i = speye(n-1);
C = [z i; -1 z'];

diff = norm(C*A - A*C, 'fro');
