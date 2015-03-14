## TxSSA ##

TxSSA is an acronym for [Tech-X Corporation](http://www.txcorp.com) Sparse Spectral Approximation. It is a library that has interfaces in C++, C, and MATLAB.  It is an implementation of a matrix sparsification algorithm that takes a general real or complex matrix as input and produces a sparse output matrix of the same size.  The non-zero entries are chosen to minimize changes to the singular values and singular vectors corresponding to the near null-space.  The output matrix is constrained to preserve left and right null-spaces exactly.  The sparsity pattern of the output matrix is automatically determined or can be given as input.

If the input matrix belongs  to a common matrix subspace, the generated sparse matrix belongs to the same subspace automatically. This holds for the subspaces of Hermitian, complex-symmetric, Hamiltonian, circulant, centrosymmetric, and persymmetric matrices, and for each of the skew counterparts.

## Authors ##

  * [Chetan Jhurani](http://www.ices.utexas.edu/~chetan/) ([chetan.jhurani@gmail.com](mailto:chetan.jhurani@gmail.com))
  * Travis M. Austin ([austin@txcorp.com](mailto:austin@txcorp.com))

[Tech-X Corporation](http://www.txcorp.com), 5621 Arapahoe Ave, Boulder, CO 80303

## A Visual ##

Darker pixels are matrix values with large magnitudes and vice versa.<br>
<img src='http://txssa.googlecode.com/svn/trunk/figures/transform.png' />

<h2>Downloads</h2>

See the source <a href='http://code.google.com/p/txssa/source/checkout'>checkout</a> page<br>
or the <a href='http://code.google.com/p/txssa/downloads/list'>download</a> page for source and binary archives.<br>
<br>
<h2>Dependencies</h2>

BLAS, LAPACK, and C and C++ compilers (C89 and C++98). CMake can be used for<br>
building the library on any supported platform.  Visual Studio projects are also provided for building<br>
on Windows.<br>
<br>
<h2>Acknowledgments</h2>

The development of this library was supported by the US Department of Energy SBIR Grant DE-FG02-08ER85154.<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<br>
<hr />