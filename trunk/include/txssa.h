/*
TxSSA: Tech-X Sparse Spectral Approximation
Copyright (C) 2012 Tech-X Corporation, 5621 Arapahoe Ave, Boulder CO 80303

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Tech-X Corporation nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
Authors:

1. Chetan Jhurani (chetan.jhurani@gmail.com, jhurani@txcorp.com)
     For more information and relevant publications, visit
     http://www.ices.utexas.edu/~chetan/

2. Travis M. Austin (austin@txcorp.com)

Contact address:

Tech-X Corporation
5621 Arapahoe Ave
Boulder, CO 80303
http://www.txcorp.com

*/


#ifndef TXSSA_H
#define TXSSA_H

/* TxSSA (Tech-X Corporation Sparse Spectral Approximation): */
/* API to compute a matrix structure preserving sparse spectral approximation */

/* This header can be included in both C and C++ translation units. */

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------- */

/*

C API Naming convention:  Function names, for example "ssa_d_pat", has following
components.

1. [ssa]
        sparse spectral approximation, characterizing the algorithm 

2. [d|s|z|c]
        d = double-real
        s = single-real
        z = double-complex
        c = single-complex

3. [pat|lpn|ids]
        pat = sparsity pattern is given by user
        lpn = L_p norm based sparse matrix to be computed by the library
        ids = L_p norm based sparsity pattern to be computed by the library

C++ API does not need [d|s|z|c], and thus the name is like ssa_pat.

*/

/* -------------------------------------------------------------------------- */

/*

Notes.
- All APIs that return an int return 0 on success and non-zero on failure.
- Use ssa_error_* APIs to get an error stack.
- Input matrix is always column-oriented, as implied by the argument names
  col_values and col_leading_dim.
- Output pattern and values are row-oriented.
- max_num_bins equal to zero means perform no binning.
- For impose_null_spaces, 0 => don't impose and 1 => impose (both left and
  right)
- sparsity_norm_p is norm used for computing pattern in [0,inf].
- sparsity_ratio is ratio of Lp norm to preserve for pattern (say 0.7-0.8).
- "reserved" pointer is used by the library.  Should not be modified or used
  by the user.
- If the library computes the output pattern and values, call deallocate when
  done.
- We don't use any language or library defined type for complex data in the
  API.  User needs to perform some casting on a real type pointer which points
  to complex data stored as (real, imag) pairs.
- If the user gives the pattern and specifies a specific matrix type, then
  the pattern must be appropriate for that matrix type.  For example, for
  Hermitian or skew-Hermitian or complex-symmetric input matrix, the input
  pattern must be symmetric.  Otherwise, because of the optimizations used
  in the code in assuming such a pattern, the result is not mathematically
  well-defined.
- One can of course not specify the matrix type and still get back the sparse
  approximation.  Every effort is made so that the sparse approximation
  retains the special structural properties if the matrix is special but
  is given as general.
- Unlike BLAS or LAPACK, if the input matrix has special structure (like
  symmetry), one must give the full matrix and not just one half of it.

*/

/* -------------------------------------------------------------------------- */

#ifdef TXSSA_API
#error TXSSA_API must not be defined before this file is included.
#endif

#if defined(_WINDOWS) || defined(_MSC_VER)

#if defined(TXSSA_DLL_EXPORTS)
#define TXSSA_API __declspec(dllexport)
#elif defined(TXSSA_DLL_IMPORTS)
#define TXSSA_API __declspec(dllimport)
#else
#define TXSSA_API
#endif

#else

#define TXSSA_API

#endif

/* -------------------------------------------------------------------------- */

/* Version numbers will be incremented with each publicly available change to
   the library interface and/or implementation.
*/

/* Header version (intf = interface). */
const unsigned int ssa_intf_version = 1;

/* Source version (impl = implementation). */
TXSSA_API extern const unsigned int ssa_impl_version;

/* -------------------------------------------------------------------------- */

/* CSR data for the case when the pattern is also computed by the library.    */
/* CSR = Compressed Sparse Row */

struct TXSSA_API ssa_d_csr
{
    int*        row_offsets;
    int*        column_ids;
    double*     values;
    const void* reserved;
};

struct TXSSA_API ssa_s_csr
{
    int*        row_offsets;
    int*        column_ids;
    float*      values;
    const void* reserved;
};

/*
(real, imag) pairs, cast to a suitable complex scalar supported in the
 calling language.
*/
struct TXSSA_API ssa_z_csr
{
    int*        row_offsets;
    int*        column_ids;
    double*     values;
    const void* reserved;
};

/*
(real, imag) pairs, cast to a suitable complex scalar supported in the
 calling language.
*/
struct TXSSA_API ssa_c_csr
{
    int*        row_offsets;
    int*        column_ids;
    float*      values;
    const void* reserved;
};

/* -------------------------------------------------------------------------- */

/* Call these to deallocate when done with the data. */

TXSSA_API void ssa_d_csr_deallocate(struct ssa_d_csr* matrix);
TXSSA_API void ssa_s_csr_deallocate(struct ssa_s_csr* matrix);
TXSSA_API void ssa_z_csr_deallocate(struct ssa_z_csr* matrix);
TXSSA_API void ssa_c_csr_deallocate(struct ssa_c_csr* matrix);

/* -------------------------------------------------------------------------- */

/* We use "Hermitian" instead of "symmetric" for real matrices too.           */
/* That is what makes sense mathematically and physically, but unfortunately
   is not used in BLAS and LAPACK. */

/* If type is unknown, use ssa_matrix_type_general.                           */

enum ssa_matrix_type
{
    /* Do not change the order. */

    ssa_matrix_type_undefined = -1,
    ssa_matrix_type_general,
    ssa_matrix_type_hermitian_pos_def,
    ssa_matrix_type_hermitian_pos_semi_def,
    ssa_matrix_type_hermitian,
    ssa_matrix_type_skew_hermitian,
    ssa_matrix_type_complex_symmetric,
    ssa_matrix_type_num_types
};

/* -------------------------------------------------------------------------- */

/* Double precision real APIs. */

/* User-given pattern */
TXSSA_API int ssa_d_pat(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    double*              out_row_values);

/* User-given parameters for computing L_p norm based matrix. */
TXSSA_API int ssa_d_lpn(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_d_csr*    out_matrix);

/* Compute the row-offsets and column-ids only, not the values. */
/* The space for values is allocated with the appropriate size. */
TXSSA_API int ssa_d_ids(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_d_csr*    out_matrix);

/* -------------------------------------------------------------------------- */

/* Single precision real APIs. */

/* User-given pattern */
TXSSA_API int ssa_s_pat(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    float*               out_row_values);

/* User-given parameters for computing L_p norm based matrix. */
TXSSA_API int ssa_s_lpn(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_s_csr*    out_matrix);

/* Compute the row-offsets and column-ids only, not the values. */
/* The space for values is allocated with the appropriate size. */
TXSSA_API int ssa_s_ids(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_s_csr*    out_matrix);

/* -------------------------------------------------------------------------- */

/* Double precision complex APIs. */

/* User-given pattern */
TXSSA_API int ssa_z_pat(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    double*              out_row_values);

/* User-given parameters for computing L_p norm based matrix. */
TXSSA_API int ssa_z_lpn(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_z_csr*    out_matrix);

/* Compute the row-offsets and column-ids only, not the values. */
/* The space for values is allocated with the appropriate size. */
TXSSA_API int ssa_z_ids(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_z_csr*    out_matrix);

/* -------------------------------------------------------------------------- */

/* Single precision complex APIs. */

/* User-given pattern */
TXSSA_API int ssa_c_pat(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    float*               out_row_values);

/* User-given parameters for computing L_p norm based matrix. */
TXSSA_API int ssa_c_lpn(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_c_csr*    out_matrix);

/* Compute the row-offsets and column-ids only, not the values. */
/* The space for values is allocated with the appropriate size. */
TXSSA_API int ssa_c_ids(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_c_csr*    out_matrix);

/* -------------------------------------------------------------------------- */

/* Error API.  Provides pointers to C strings corresponding to errors.        */
/* Each non ssa_error_* API clears the stack when called.                     */

/* -1 return value means something went wrong in error API */

/* Returns number of entries in error stack. */
TXSSA_API int ssa_error_size();

/* Get pointer to a C string for ith error.  0 <= i < ssa_error_size(). */
/* Low i means low-level (child) function and vice versa. */
TXSSA_API int ssa_error_string(int i, const char** ptr_to_error_string);

/* Clears the error stack. */
TXSSA_API int ssa_error_clear();

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus
} // extern "C"
#endif

/* -------------------------------------------------------------------------- */

#if defined(__cplusplus) && !defined(TXSSA_NO_CPP_API)

#include <complex>

// C++ APIs

// Notes on template types:
// - index_type can be any signed or unsigned, int or long long type.
// - offset_type can be any integral of same sign as index_type and not smaller
//   than index_type.
// - value_type can be float or double.
// - scalar_type can be float or double.  It is used to signifiy the precision
//   for APIs for complex data type.

// To get the pattern and values, user can default-construct a ssa_csr
// object and it will be cleaned up when the destructor is called.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
class TXSSA_API ssa_csr
{
public:

    ssa_csr();
    ~ssa_csr();

    offset_type* row_offsets;
    index_type*  column_ids;
    value_type*  values;
    const void*  reserved;

private:

    // No need for user to mess with these.
    ssa_csr(const ssa_csr&);
    ssa_csr& operator=(const ssa_csr&);
};

// -----------------------------------------------------------------------------

// real and complex API.  value_type can be real or complex.

// User-given pattern.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
TXSSA_API int ssa_pat(
    index_type            num_rows,
    index_type            num_cols,
    const value_type*     col_values,
    index_type            col_leading_dim,
    const offset_type*    row_offsets,
    const index_type*     column_ids,
    offset_type           max_num_bins,
    bool                  impose_null_spaces,
    enum ssa_matrix_type  matrix_type,
    value_type*           out_row_values);

// -----------------------------------------------------------------------------

// real API. value_type can be real only.

// User-given parameters for computing L_p norm based matrix.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
TXSSA_API int ssa_lpn(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    value_type           sparsity_ratio,
    value_type           sparsity_norm_p,
    offset_type          max_num_bins,
    bool                 impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix);

// Compute the row-offsets and column-ids only, not the values.
// The space for values is allocated with the appropriate size.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
TXSSA_API int ssa_ids(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    value_type           sparsity_ratio,
    value_type           sparsity_norm_p,
    index_type           min_num_nnz_per_row,
    index_type           min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix);

// -----------------------------------------------------------------------------

// std::complex<scalar_type> API.  scalar_type can be real only.

// User-given parameters for computing L_p norm based matrix.
template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
TXSSA_API int ssa_lpn(
    index_type                       num_rows,
    index_type                       num_cols,
    const std::complex<scalar_type>* col_values,
    index_type                       col_leading_dim,
    scalar_type                      sparsity_ratio,
    scalar_type                      sparsity_norm_p,
    offset_type                      max_num_bins,
    bool                             impose_null_spaces,
    enum ssa_matrix_type             matrix_type,
    ssa_csr<index_type, offset_type, std::complex<scalar_type> >& out_matrix);

// Compute the row-offsets and column-ids only, not the values.
// The space for values is allocated with the appropriate size.
template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
TXSSA_API int ssa_ids(
    index_type           num_rows,
    index_type           num_cols,
    const std::complex<scalar_type>* col_values,
    index_type           col_leading_dim,
    scalar_type          sparsity_ratio,
    scalar_type          sparsity_norm_p,
    index_type           min_num_nnz_per_row,
    index_type           min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, std::complex<scalar_type> >& out_matrix);

// -----------------------------------------------------------------------------

#endif /* __cplusplus */

/* -------------------------------------------------------------------------- */

#endif /* TXSSA_H */
