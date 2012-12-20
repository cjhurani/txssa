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


#ifndef DENSE_MATRIX_PERMUTE_H
#define DENSE_MATRIX_PERMUTE_H

// -----------------------------------------------------------------------------

#include "blas_wrap/dense_vector_utils.h"
#include "math/vector_utils.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_permute_rows(
    index_type  num_rows,
    index_type  num_cols,
    const value_type* A_col_values,
    index_type  A_col_leading_dim,
    const index_type* pivots,   // size num_rows
    index_type  pivots_base,    // 1 if coming from LAPACK.
    value_type* row_permuted_A, // num_rows x num_cols
    index_type  row_permuted_A_leading_dim)
{
    bool success =
        A_col_values &&
        pivots &&
        row_permuted_A &&
        num_rows <= A_col_leading_dim &&
        num_rows <= row_permuted_A_leading_dim;

    if(success)
    {
        // pivots[j] = i + pivots_base iff P(i,j) = delta_ij

        const index_type copy_r_n    = num_cols;
        const index_type copy_r_incx = A_col_leading_dim;
        const index_type copy_r_incy = row_permuted_A_leading_dim;

        for(index_type step = 0; (step < num_rows) && success; ++step)
        {
            const value_type* copy_r_x = A_col_values + step;
            value_type* copy_r_y = row_permuted_A + pivots[step] - pivots_base;

            success = dense_vector_utils_copy(copy_r_n, copy_r_x, copy_r_incx, copy_r_y, copy_r_incy);
            assert(success);
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_permute_rows: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_permute_rows_blocked(
    index_type  num_rows,
    index_type  num_cols,
    const value_type* A_col_values,
    index_type  A_col_leading_dim,
    const index_type* pivots,   // size num_rows
    index_type  pivots_base,    // 1 if coming from LAPACK.
    value_type* row_permuted_A, // num_rows x num_cols
    index_type  row_permuted_A_leading_dim,
    index_type  block_size) // in [1, num_cols]
{
    bool success =
        0 < block_size &&
        block_size <= num_cols &&
        A_col_values &&
        pivots &&
        num_rows <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const value_type* A_col_it = A_col_values;
        value_type* permuted_A_it = row_permuted_A;

        for(index_type col = 0; (col < num_cols) && success; col += block_size)
        {
            const index_type this_block_size =
                (num_cols < col + block_size) ? num_cols - col : block_size;

            success = dense_matrix_permute_rows(
                num_rows,
                this_block_size,
                A_col_it,
                A_col_leading_dim,
                pivots,
                pivots_base,
                permuted_A_it,
                row_permuted_A_leading_dim);

            assert(success);

            A_col_it += std::size_t(this_block_size) * std::size_t(A_col_leading_dim);
            permuted_A_it += std::size_t(this_block_size) * std::size_t(row_permuted_A_leading_dim);
        }
    }

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_permute_rows_blocked: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_permute_cols(
    index_type  num_rows,
    index_type  num_cols,
    const value_type* A_col_values,
    index_type  A_col_leading_dim,
    const index_type* pivots,   // size num_cols
    index_type  pivots_base,    // 1 if coming from LAPACK.
    value_type* col_permuted_A, // num_rows x num_cols
    index_type  col_permuted_A_leading_dim)
{
    bool success =
        A_col_values &&
        pivots &&
        col_permuted_A &&
        num_rows <= A_col_leading_dim &&
        num_rows <= col_permuted_A_leading_dim;

    assert(success);

    if(success)
    {
        // pivots[j] = i + pivots_base iff P(i,j) = delta_ij

        const index_type copy_c_n    = num_rows;
        const index_type copy_c_incx = 1;
        const index_type copy_c_incy = 1;

        for(index_type step = 0; (step < num_cols) && success; ++step)
        {
            const value_type* copy_c_x = A_col_values + step * A_col_leading_dim;
            value_type* copy_c_y       = col_permuted_A + (pivots[step] - pivots_base) * col_permuted_A_leading_dim;

            success = dense_vector_utils_copy(copy_c_n, copy_c_x, copy_c_incx, copy_c_y, copy_c_incy);
            assert(success);
        }
    }

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_permute_cols: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Overwrite a square A by its row and column permutation.  Use given
// block size for row permutations.

template<typename index_type, typename value_type>
bool dense_matrix_permute_sym_blocked(
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    const index_type* pivots,  // size matrix_size
    index_type  pivots_base,   // 1 if coming from LAPACK.
    value_type* work,          // matrix_size^2
    index_type col_block_size)
{
    // pivots[j] = i + pivots_base iff P(i,j) = delta_ij

    bool success =
        A_col_values &&
        pivots &&
        work &&
        matrix_size <= A_col_leading_dim
        &&
        dense_matrix_permute_rows_blocked(
            matrix_size, matrix_size,
            A_col_values, A_col_leading_dim,
            pivots, pivots_base,
            work, matrix_size, col_block_size)
        &&
        dense_matrix_permute_cols(
            matrix_size, matrix_size,
            work, matrix_size,
            pivots, pivots_base,
            A_col_values, A_col_leading_dim);

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_permute_sym_blocked: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Overwrite a square A by its row and column permutation.

template<typename index_type, typename value_type>
bool dense_matrix_permute_sym(
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    const index_type* pivots,  // size matrix_size
    index_type  pivots_base,   // 1 if coming from LAPACK.
    value_type* work)          // matrix_size^2
{
    bool success = dense_matrix_permute_sym_blocked(
        matrix_size,
        A_col_values,
        A_col_leading_dim,
        pivots,
        pivots_base,
        work,
        matrix_size);

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_permute_sym: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_PERMUTE_H
