function diff = diff_hamiltonian(A)

% diff = diff_hamiltonian(A)
%
% Checks how far is a matrix A different from being Hamiltonian.
% A Frobenius norm based difference is reported.

if(ndims(A) ~= 2 || size(A,1) ~= size(A,2) || mod(size(A,1),2) ~= 0)
    assert(false, 'diff_hamiltonian: Input A is not 2-D or not square or not even-sized');
end

K = spkay(size(A,1));

diff = norm(K*A + A'*K, 'fro');

