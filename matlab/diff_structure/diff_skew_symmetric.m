function diff = diff_skew_symmetric(A)

% diff = diff_skew_symmetric(A)
%
% Checks how far is a real matrix A different from being skew_symmetric.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2) || ~isreal(A))
    assert(false, 'diff_skew_symmetric: Input A is not 2-D or not square or not real');
end

diff = norm(A + A','fro');

