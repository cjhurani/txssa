Repository migrated from https://code.google.com/p/txssa

TxSSA stands for Tech-X Corporation Sparse Spectral Approximation.

The library TxSSA, which has interfaces in C++, C, and MATLAB, is an
implementation of a newly designed matrix sparsification algorithm
that takes a general real or complex matrix as input and produces a
sparse output matrix of the same size. The non-zero entries are chosen
to minimize changes to the singular values and singular vectors
corresponding to the near null-space. The output matrix is constrained
to preserve left and right null-spaces exactly. The sparsity pattern
of the output matrix is automatically determined or can be given as
input.

If the input matrix belongs to a common matrix subspace, the generated
sparse matrix belongs to the same subspace automatically. This holds
for the subspaces of Hermitian, complex-symmetric, Hamiltonian,
circulant, centrosymmetric, and persymmetric matrices, and for each of
the skew counterparts.

Notes:

- The development of this library was supported by the US Department of
  Energy SBIR Grant DE-FG02-08ER85154.

- The MATLAB code is in ./matlab directory.  The library with C/C++ API can
  be built using CMake.  See install.txt for instructions.

- Version numbers are like x.y where x >= 1 is the interface version and
  y >= 0 is the implementation version.
