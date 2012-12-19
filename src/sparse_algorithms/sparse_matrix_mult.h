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


#ifndef SPARSE_MATRIX_MULT_H
#define SPARSE_MATRIX_MULT_H

// -----------------------------------------------------------------------------

#include "dense_vectors/dense_vectors_utils.h"
#include "dense_vectors/dense_vectors.h"
#include "sparse_vectors/sparse_vectors.h"
#include "math/complex_types.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Multiplication of sparse vectors with vectors.
// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_matrix_mult(
    index_type num_rows,
    index_type num_cols,
    const offset_type* offsets,
    const index_type* ids,
    const value_type* values,
    index_type x_num_cols,
    const value_type* x,         // num_cols x x_num_cols
    index_type x_col_leading_dim,
    value_type* Ax,              // num_rows x x_num_cols
    index_type Ax_col_leading_dim,
    value_type mult_factor)
{
    bool success =
        offsets &&
        ids &&
        values &&
        x &&
        Ax &&
        num_cols <= x_col_leading_dim &&
        num_rows <= Ax_col_leading_dim;

    assert(success);

    if(success)
    {
        const value_type* xk = x;
        value_type* Axk = Ax;

        for(index_type k = 0; k < x_num_cols; ++k)
        {
            for(index_type i = 0; i < num_rows; ++i)
            {
                value_type dotp = value_type();

                for(offset_type jj = offsets[i]; jj < offsets[i+1]; ++jj)
                {
                    assert(ids[jj] < num_cols);
                    dotp += values[jj] * xk[ids[jj]];
                }

                Axk[i] = mult_factor * dotp;
            }

            xk += x_col_leading_dim;
            Axk += Ax_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "sparse_matrix_mult: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_matrix_mult_trans(
    index_type num_rows,
    index_type num_cols,
    const offset_type* offsets,
    const index_type* ids,
    const value_type* values,
    index_type x_num_cols,
    const value_type* x,          // num_rows x x_num_cols
    index_type x_col_leading_dim,
    value_type* ATx,              // num_cols x x_num_cols
    index_type ATx_col_leading_dim,
    value_type mult_factor)
{
    bool success =
        offsets &&
        ids &&
        values &&
        x &&
        ATx &&
        num_rows <= x_col_leading_dim &&
        num_cols <= ATx_col_leading_dim &&
        dense_vectors_utils_fill(
            x_num_cols,
            num_cols,
            ATx,
            ATx_col_leading_dim,
            value_type(0));

    assert(success);

    if(success)
    {
        for(index_type i = 0; i < num_rows; ++i)
        {
            const value_type* x_it = x;
            value_type* ATx_it = ATx;

            for(index_type x_col = 0; x_col < x_num_cols; ++x_col)
            {
                const value_type factor_x_val = mult_factor * x_it[i];

                for(offset_type jj = offsets[i]; jj < offsets[i+1]; ++jj)
                {
                    assert(ids[jj] < num_cols);
                    ATx_it[ids[jj]] += factor_x_val * std::conj(values[jj]);
                }

                x_it += x_col_leading_dim;
                ATx_it += ATx_col_leading_dim;
            }
        }
    }

    if(!success)
        internal_api_error_set_last(
            "sparse_matrix_mult_trans: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type,
    typename value_type_x,
    typename value_type_Ax
>
bool sparse_matrix_mult(
    const sparse_vectors<index_type, offset_type, value_type>& A_row,
    const dense_vectors<index_type, value_type_x>& x_col,
          dense_vectors<index_type, value_type_Ax>& Ax_col,
    value_type mult_factor)
{
    bool success =
        A_row.max_size() == x_col.max_size() &&
        x_col.num_vecs() == Ax_col.num_vecs() &&
        A_row.num_vecs() == Ax_col.max_size() &&
        sparse_matrix_mult(
            A_row.num_vecs(),
            A_row.max_size(),
            A_row.vec_offsets(),
            A_row.vec_ids(),
            A_row.vec_values(),
            x_col.num_vecs(),
            x_col.vec_values(),
            x_col.leading_dimension(),
            Ax_col.vec_values(),
            Ax_col.leading_dimension(),
            mult_factor);

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "sparse_matrix_mult: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type,
    typename value_type_x,
    typename value_type_ATx
>
bool sparse_matrix_mult_trans(
    const sparse_vectors<index_type, offset_type, value_type>& A_row,
    const dense_vectors<index_type, value_type_x>& x_col,
          dense_vectors<index_type, value_type_ATx>& ATx_col,
    value_type mult_factor)
{
    bool success =
        A_row.max_size() == ATx_col.max_size() &&
        x_col.num_vecs() == ATx_col.num_vecs() &&
        A_row.num_vecs() == x_col.max_size() &&
        sparse_matrix_mult_trans(
            A_row.num_vecs(),
            A_row.max_size(),
            A_row.vec_offsets(),
            A_row.vec_ids(),
            A_row.vec_values(),
            x_col.num_vecs(),
            x_col.vec_values(),
            x_col.leading_dimension(),
            ATx_col.vec_values(),
            ATx_col.leading_dimension(),
            mult_factor);

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "sparse_matrix_mult_trans: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_MATRIX_MULT_H
