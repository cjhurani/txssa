function diff = diff_persymmetric(A)

% diff = diff_persymmetric(A)
%
% Checks how far is a matrix A different from being persymmetric.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2))
    assert(false, 'diff_persymmetric: Input A is not 2-D or not square');
end

J = spjay(size(A,1));

diff = norm(A*J - J*A','fro');

