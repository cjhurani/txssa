function diff = diff_hermitian(A)

% diff = diff_hermitian(A)
%
% Checks how far is a matrix A different from being hermitian.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2))
    assert(false, 'diff_hermitian: Input A is not 2-D or not square');
end

diff = norm(A - A','fro');

