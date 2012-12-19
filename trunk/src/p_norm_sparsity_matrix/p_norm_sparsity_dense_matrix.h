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


#ifndef P_NORM_SPARSITY_DENSE_MATRIX_H
#define P_NORM_SPARSITY_DENSE_MATRIX_H

// -----------------------------------------------------------------------------

#include "sparsity_union/sparse_vectors_union_w_trans.h"
#include "p_norm_sparsity_vectors/p_norm_sparsity_dense_vectors.h"
#include "dense_vectors/dense_vectors_transpose_view.h"
#include "dense_vectors/dense_vectors.h"
#include "cpp/vector_vector_id.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <sstream>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: To compute the p-norm sparsity pattern of dense matrices given in
// various forms.
// 1. Both row and col oriented values given for a general matrix.  They would
//    refer to the same numerical matrix.  The purpose is that if data in both
//    orientations is known, then it could be cache-friendly to use what is
//    given rather than make large jumps to access the transpose.
// 2. Row or col values given for a matrix whose absolute value is symmetric.
// 3. Row-oriented general dense matrix.
// 3. Col-oriented general dense matrix.
// -----------------------------------------------------------------------------

// Output ids will be sorted and unique.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_matrix(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz_per_row,
    index_type min_num_nnz_per_col,

// input data:
    index_type num_rows,
    index_type num_cols,
    const value_type* row_values,
    index_type row_leading_dim,
    const value_type* col_values,
    index_type col_leading_dim,

// output ids:
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(
        !row_values ||
        !col_values ||
        row_leading_dim < num_cols ||
        col_leading_dim < num_rows)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_matrix: Unacceptable input argument(s).";
        if(!row_values)                oss << " !row_values.";
        if(!col_values)                oss << " !col_values.";
        if(row_leading_dim < num_cols) oss << " row_leading_dim < num_cols.";
        if(col_leading_dim < num_rows) oss << " col_leading_dim < num_rows.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_row_pat, tmp_col_pat;

    bool success =
        // row-wise pattern
        p_norm_sparsity_dense_vectors(
            ratio,
            p,
            min_num_nnz_per_row,
            num_rows,
            num_cols,
            row_leading_dim,
            row_values,
            tmp_row_pat)
        &&
        // col-wise pattern
        p_norm_sparsity_dense_vectors(
            ratio,
            p,
            min_num_nnz_per_col,
            num_cols,
            num_rows,
            col_leading_dim,
            col_values,
            tmp_col_pat);

    // ids in tmp_row_pat and tmp_col_pat are not sorted necessarily.

    if(success)
    {
        const dense_vectors<index_type, const value_type> row_matrix(
            num_rows,
            num_cols,
            row_leading_dim,
            row_values);

        // union of the two patterns
        success = sparse_vectors_union_w_trans(
            vector_vector_id<index_type>(tmp_row_pat),
            row_matrix,
            vector_vector_id<index_type>(tmp_col_pat),
            row_oriented_sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_matrix: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// If matrix is actually abs_sym, this will work for both row- and column-
// oriented data.  Output ids will be sorted and unique.
// "abs_sym" means its element-wise absolute value is symmetric.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_matrix_abs_sym(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz,

// input data:
    index_type matrix_size,
    const value_type* values,
    index_type leading_dim,

// output ids:
    std::vector< std::vector<index_type> >& sparse_pat)
{
    if(
        !values ||
        leading_dim < matrix_size)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_matrix_abs_sym: Unacceptable input argument(s).";
        if(!values)                   oss << " !values.";
        if(leading_dim < matrix_size) oss << " leading_dim < matrix_size.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_pat;

    bool success =
        p_norm_sparsity_dense_vectors(
            ratio,
            p,
            min_num_nnz,
            matrix_size,
            matrix_size,
            leading_dim,
            values,
            tmp_pat);

    if(success)
    {
        const dense_vectors<index_type, const value_type> the_matrix(
            matrix_size,
            matrix_size,
            leading_dim,
            values);

        success = sparse_vectors_union_w_self_trans(
            vector_vector_id<index_type>(tmp_pat),
            the_matrix,
            sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_matrix_abs_sym: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Output ids will be sorted and unique.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_matrix_row_oriented(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz_per_row,
    index_type min_num_nnz_per_col,

// input data:
    index_type num_rows,
    index_type num_cols,
    const value_type* row_values, // values in rows
    index_type row_leading_dim,

// output ids:
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(
        !row_values ||
        row_leading_dim < num_cols)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_matrix_row_oriented: Unacceptable input argument(s).";
        if(!row_values)                oss << " !row_values.";
        if(row_leading_dim < num_cols) oss << " row_leading_dim < num_cols.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_row_pat, tmp_col_pat;

    bool success =
        // row-wise pattern
        p_norm_sparsity_dense_vectors(
            ratio,
            p,
            min_num_nnz_per_row,
            num_rows,
            num_cols,
            row_leading_dim,
            row_values,
            tmp_row_pat)
        &&
        // col-wise pattern
        p_norm_sparsity_dense_vectors_transpose_view(
            ratio,
            p,
            min_num_nnz_per_col,
            num_rows,
            num_cols,
            row_leading_dim,
            row_values,
            tmp_col_pat);

    // ids in tmp_row_pat and tmp_col_pat are not sorted necessarily.

    if(success)
    {
        const dense_vectors<index_type, const value_type> row_matrix(
            num_rows,
            num_cols,
            row_leading_dim,
            row_values);

        // union of the two patterns
        success = sparse_vectors_union_w_trans(
            vector_vector_id<index_type>(tmp_row_pat),
            row_matrix,
            vector_vector_id<index_type>(tmp_col_pat),
            row_oriented_sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_matrix_row_oriented: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Output ids will be sorted and unique.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_matrix_col_oriented(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz_per_row,
    index_type min_num_nnz_per_col,

// input data:
    index_type num_rows,
    index_type num_cols,
    const value_type* col_values, // values in cols
    index_type col_leading_dim,

// output ids:
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(
        !col_values ||
        col_leading_dim < num_rows)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_matrix_col_oriented: Unacceptable input argument(s).";
        if(!col_values)                oss << " !col_values.";
        if(col_leading_dim < num_cols) oss << " col_leading_dim < num_cols.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_row_pat, tmp_col_pat;

    bool success =
        // row-wise pattern
        p_norm_sparsity_dense_vectors_transpose_view(
            ratio,
            p,
            min_num_nnz_per_row,
            num_cols,
            num_rows,
            col_leading_dim,
            col_values,
            tmp_row_pat)
        &&
        // col-wise pattern
        p_norm_sparsity_dense_vectors(
            ratio,
            p,
            min_num_nnz_per_col,
            num_cols,
            num_rows,
            col_leading_dim,
            col_values,
            tmp_col_pat);

    // ids in tmp_row_pat and tmp_col_pat are not sorted necessarily.

    if(success)
    {
        const dense_vectors<index_type, const value_type> col_matrix(
            num_cols,
            num_rows,
            col_leading_dim,
            col_values);

        const dense_vectors_transpose_view<index_type, const value_type> row_matrix(
            col_matrix);

        // union of the two patterns
        success = sparse_vectors_union_w_trans(
            vector_vector_id<index_type>(tmp_row_pat),
            row_matrix,
            vector_vector_id<index_type>(tmp_col_pat),
            row_oriented_sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_matrix_col_oriented: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_DENSE_MATRIX_H
