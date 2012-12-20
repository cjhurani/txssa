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


#ifndef P_NORM_SPARSITY_SPARSE_MATRIX_H
#define P_NORM_SPARSITY_SPARSE_MATRIX_H

// -----------------------------------------------------------------------------

#include "sparsity_union/sparse_vectors_union_w_trans.h"
#include "p_norm_sparsity_vectors/p_norm_sparsity_sparse_vectors.h"
#include "sparse_vectors/sparse_vectors.h"
#include "math/precision_traits.h"
#include "cpp/vector_vector_id.h"
#include "internal_api_error/internal_api_error.h"
#include <sstream>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: To compute the p-norm sparsity pattern of sparse matrices given in
// various forms.
// 1. Both row and col oriented values given for a general matrix.  They would
//    refer to the same numerical matrix.  The purpose is that if data in both
//    orientations is known, then it would be easier to use what is given rather
//    than create a transpose.
// 2. Row or col values given for a matrix whose absolute value is symmetric.
// 3. Row-oriented general sparse matrix, for which we'll have to create a
//    transpose first.
// -----------------------------------------------------------------------------

// Output ids will be sorted, unique, and relative to row ids.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_norm_sparsity_sparse_matrix(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz_per_row,
    index_type min_num_nnz_per_col,

// input data:
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets, // offsets in rows
    const index_type* row_ids,      // ids in rows, must be sorted.
    const value_type* row_values,   // values in rows
    const offset_type* col_offsets, // offsets in columns
    const index_type* col_ids,      // ids in columns, need not be sorted.
    const value_type* col_values,   // values in columns

// output ids:
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(
        !row_offsets || !row_ids || !row_values ||
        !col_offsets || !col_ids || !col_values)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_sparse_matrix: Unacceptable input argument(s).";
        if(!row_offsets) oss << " !row_offsets.";
        if(!row_ids)     oss << " !row_ids.";
        if(!row_values)  oss << " !row_values.";
        if(!col_offsets) oss << " !col_offsets.";
        if(!col_ids)     oss << " !col_ids.";
        if(!col_values)  oss << " !col_values.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_row_pat, tmp_col_pat;

    bool success =
        // row-wise pattern
        p_norm_sparsity_sparse_vectors(
            ratio,
            p,
            min_num_nnz_per_row,
            num_rows,
            num_cols,
            row_offsets,
            row_values,
            tmp_row_pat)
        &&
        // col-wise pattern
        p_norm_sparsity_sparse_vectors(
            ratio,
            p,
            min_num_nnz_per_col,
            num_cols,
            num_rows,
            col_offsets,
            col_values,
            tmp_col_pat);

    if(success)
    {
        sparse_vectors_ids<const index_type, const offset_type>
            row_id_vecs(num_rows, num_cols, row_offsets, row_ids);

        // union of the two patterns
        success = sparse_vectors_union_w_trans(
            vector_vector_id<index_type>(tmp_row_pat),
            row_id_vecs,
            vector_vector_id<index_type>(tmp_col_pat),
            row_oriented_sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_sparse_matrix: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// If matrix is actually abs_sym, this will work for both row- and column-
// oriented data.  Output ids will be sorted, unique, and relative to row ids.
// "abs_sym" means its element-wise absolute value is symmetric.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_norm_sparsity_sparse_matrix_abs_sym(
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz,
    index_type matrix_size,
    const offset_type* offsets, // offsets in rows/cols
    const index_type* ids,      // ids in rows/cols, must be sorted
    const value_type* values,   // values in rows/cols
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(!offsets || !ids || !values)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_sparse_matrix_abs_sym: Unacceptable input argument(s).";
        if(!offsets) oss << " !offsets.";
        if(!ids)     oss << " !ids.";
        if(!values)  oss << " !values.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    std::vector< std::vector<index_type> > tmp_pat;

    bool success =
        p_norm_sparsity_sparse_vectors(
            ratio,
            p,
            min_num_nnz,
            matrix_size,
            matrix_size,
            offsets,
            values,
            tmp_pat);

    if(success)
    {
        sparse_vectors_ids<const index_type, const offset_type>
            id_vecs(matrix_size, matrix_size, offsets, ids);

        success = sparse_vectors_union_w_self_trans(
            vector_vector_id<index_type>(tmp_pat),
            id_vecs,
            row_oriented_sparse_pat);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_sparse_matrix_abs_sym: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Output ids will be sorted, unique, and relative to row ids.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_norm_sparsity_sparse_matrix_row_oriented(
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz_per_row,
    index_type min_num_nnz_per_col,
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets, // offsets in rows
    const index_type* row_ids,      // ids in rows, must be sorted.
    const value_type* row_values,   // values in rows
    std::vector< std::vector<index_type> >& row_oriented_sparse_pat)
{
    if(!row_offsets || !row_ids || !row_values)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_sparse_matrix_row_oriented: Unacceptable input argument(s).";
        if(!row_offsets) oss << " !row_offsets.";
        if(!row_ids)     oss << " !row_ids.";
        if(!row_values)  oss << " !row_values.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    const sparse_vectors
        <const index_type, const offset_type, const value_type>
        row_oriented(num_rows, num_cols, row_offsets, row_ids, row_values);

    sparse_vectors<index_type, offset_type, value_type> col_oriented;

    bool success =
        row_oriented.get_transpose(col_oriented)
        &&
        p_norm_sparsity_sparse_matrix(
                ratio,
                p,
                min_num_nnz_per_row,
                min_num_nnz_per_col,
                num_rows,
                num_cols,
                row_offsets,
                row_ids,
                row_values,
                col_oriented.vec_offsets(),
                col_oriented.vec_ids(),
                col_oriented.vec_values(),
                row_oriented_sparse_pat);

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_sparse_matrix_row_oriented: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_SPARSE_MATRIX_H
