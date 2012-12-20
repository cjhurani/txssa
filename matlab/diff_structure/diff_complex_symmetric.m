function diff = diff_complex_symmetric(A)

% diff = diff_complex_symmetric(A)
%
% Checks how far is a matrix A different from being complex_symmetric.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2))
    assert(false, 'diff_complex_symmetric: Input A is not 2-D or not square or is real');
end

diff = norm(A - A.','fro');

