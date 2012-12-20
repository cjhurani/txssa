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


#ifndef P_NORMALIZE_SPARSE_MATRIX_H
#define P_NORMALIZE_SPARSE_MATRIX_H

// -----------------------------------------------------------------------------

#include "p_normalize_vectors_and_trans.h"
#include "sparse_vectors/sparse_vectors.h"
#include "math/precision_traits.h"
#include <cstddef>

// -----------------------------------------------------------------------------
// Objective: Iteratively normalize a sparse matrix so that each row/col has
// p-norm equal to 1.
// -----------------------------------------------------------------------------

// The matrix is row-oriented, but normalization happens to both rows/cols.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_normalize_sparse_matrix_row_oriented(
// algorithmic options:
    typename precision_traits<value_type>::scalar p,
    typename precision_traits<value_type>::scalar tolerance,
    std::size_t max_iters,

// input data:
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets,
    const index_type* row_ids,

// input/output data:
    value_type* row_values,

// output data:
    typename precision_traits<value_type>::scalar* left_diag,
    typename precision_traits<value_type>::scalar* right_diag,
    std::size_t* iters)
{
    typedef sparse_vectors<const index_type, const offset_type, value_type>
        sparse_matrix_type;

    sparse_matrix_type row_oriented_sparse_matrix(
        num_rows,
        num_cols,
        row_offsets,
        row_ids,
        row_values);

    return p_normalize_vectors_and_trans
        <index_type, value_type, sparse_matrix_type>(
            p,
            tolerance,
            max_iters,
            row_oriented_sparse_matrix,
            left_diag,
            right_diag,
            iters);
}

// -----------------------------------------------------------------------------

// Matrix orientation is immaterial if matrix is abs_sym.
// "abs_sym" means the matrix's element-wise absolute value is symmetric.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_normalize_sparse_matrix_abs_sym(
// algorithmic options:
    typename precision_traits<value_type>::scalar p,
    typename precision_traits<value_type>::scalar tolerance,
    std::size_t max_iters,

// input data:
    index_type matrix_size,
    const offset_type* offsets,
    const index_type* ids,

// input/output data:
    value_type* values,

// output data:
    typename precision_traits<value_type>::scalar* diag,
    std::size_t* iters)
{
    typedef sparse_vectors<const index_type, const offset_type, value_type>
        sparse_matrix_type;

    sparse_matrix_type sparse_matrix(
        matrix_size,
        matrix_size,
        offsets,
        ids,
        values);

    return p_normalize_vectors_and_trans_abs_sym
        <index_type, value_type, sparse_matrix_type>(
            p,
            tolerance,
            max_iters,
            sparse_matrix,
            diag,
            iters);
}

// -----------------------------------------------------------------------------

#endif // P_NORMALIZE_SPARSE_MATRIX_H
