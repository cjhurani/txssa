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


#ifndef SPARSE_SPECTRAL_BINNING_H
#define SPARSE_SPECTRAL_BINNING_H

// -----------------------------------------------------------------------------

#include "matrix_binning/matrix_binning.h"
#include "matrix_binning/split_pattern_to_bins.h"
#include "sparse_vectors/sparse_vectors.h"
#include "dense_vectors/dense_vectors.h"
#include "dense_vectors/dense_vectors_transpose_view.h"
#include "internal_api_error/internal_api_error.h"
#include <limits>
#include <stdexcept>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute the split row and column bins from given pattern, matrix,
// and binning parameter.
// -----------------------------------------------------------------------------

// For real matrices.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_spectral_binning_row(
    index_type num_rows,
    index_type num_cols,
    const value_type* A_col_values,    // num_rows x num_cols
    index_type A_col_leading_dim,
    const offset_type* row_offsets,    // Size num_rows + 1
    const index_type* column_ids,      // Size row_offsets[num_rows]
    offset_type max_num_bins,
    offset_type& actual_num_bins,
    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >& row_split_pattern,
    offset_type* row_bin_ids) // Size row_offsets[num_rows]
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim &&
        row_offsets &&
        column_ids &&
        (row_bin_ids || !num_rows);

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_row: Unacceptable input argument(s).");

        return false;
    }

    std::vector<offset_type> bin_work_array;

    try
    {
        bin_work_array.resize(max_num_bins);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("sparse_spectral_binning_row: Exception. ") + exc.what()));

        return false;
    }

    const dense_vectors<index_type, const value_type> col_matrix(
        num_cols,
        num_rows,
        A_col_leading_dim,
        A_col_values);

    const dense_vectors_transpose_view<index_type, const value_type> row_matrix(
        col_matrix);

    sparse_vectors_ids<const index_type, const offset_type> in_row_pattern(
        num_rows, num_cols, row_offsets, column_ids);

    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >
        tmp_row_split_pattern;

    success =
        matrix_binning(
            in_row_pattern,
            row_matrix,
            max_num_bins,
            row_bin_ids,
            actual_num_bins,
            bin_work_array.size() ? &bin_work_array.front() : 0)
        &&
        tmp_row_split_pattern.allocate(
            actual_num_bins, offset_type(1))
        &&
        split_pattern_to_bins(
            num_rows, num_cols,
            row_offsets, column_ids,
            row_bin_ids,
            tmp_row_split_pattern.num_vecs(),
            tmp_row_split_pattern.vec_values());

    if(success)
    {
        row_split_pattern.swap(tmp_row_split_pattern);
    }
    else
    {
        actual_num_bins = std::numeric_limits<offset_type>::max();
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_row: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

// For complex matrices.
template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
bool sparse_spectral_binning_row(
    index_type num_rows,
    index_type num_cols,
    const std::complex<scalar_type>* A_col_values,    // num_rows x num_cols
    index_type A_col_leading_dim,
    const offset_type* row_offsets,    // Size num_rows + 1
    const index_type* column_ids,      // Size row_offsets[num_rows]
    offset_type max_num_bins,
    offset_type& real_actual_num_bins,
    offset_type& imag_actual_num_bins,
    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >& real_row_split_pattern,
    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >& imag_row_split_pattern,
    offset_type* real_row_bin_ids, // Size row_offsets[num_rows]
    offset_type* imag_row_bin_ids) // Size row_offsets[num_rows]
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim &&
        row_offsets &&
        column_ids &&
        real_row_bin_ids &&
        imag_row_bin_ids;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_row: Unacceptable input argument(s).");

        return false;
    }

    std::vector<offset_type> bin_work_array;

    try
    {
        bin_work_array.resize(max_num_bins);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("sparse_spectral_binning_row: Exception. ") + exc.what()));

        return false;
    }

    const dense_vectors<index_type, const std::complex<scalar_type> > col_matrix(
        num_cols,
        num_rows,
        A_col_leading_dim,
        A_col_values);

    const dense_vectors_transpose_view<index_type, const std::complex<scalar_type> > row_matrix(
        col_matrix);

    sparse_vectors_ids<const index_type, const offset_type> in_row_pattern(
        num_rows, num_cols, row_offsets, column_ids);

    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >
        tmp_real_row_split_pattern,
        tmp_imag_row_split_pattern;

    success =
        matrix_binning(
            in_row_pattern,
            row_matrix,
            max_num_bins,
            real_row_bin_ids, imag_row_bin_ids,
            real_actual_num_bins, imag_actual_num_bins,
            bin_work_array.size() ? &bin_work_array.front() : 0)
        &&
        tmp_real_row_split_pattern.allocate(
            real_actual_num_bins, offset_type(1))
        &&
        tmp_imag_row_split_pattern.allocate(
            imag_actual_num_bins, offset_type(1))
        &&
        split_pattern_to_bins(
            num_rows, num_cols,
            row_offsets,
            column_ids,
            real_row_bin_ids,
            tmp_real_row_split_pattern.num_vecs(),
            tmp_real_row_split_pattern.vec_values())
        &&
        split_pattern_to_bins(
            num_rows, num_cols,
            row_offsets,
            column_ids,
            imag_row_bin_ids,
            tmp_imag_row_split_pattern.num_vecs(),
            tmp_imag_row_split_pattern.vec_values());

    if(success)
    {
        real_row_split_pattern.swap(tmp_real_row_split_pattern);
        imag_row_split_pattern.swap(tmp_imag_row_split_pattern);
    }
    else
    {
        real_actual_num_bins = std::numeric_limits<offset_type>::max();
        imag_actual_num_bins = std::numeric_limits<offset_type>::max();
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_row: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename offset_type>
bool sparse_spectral_binning_to_col(
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets,    // Size num_rows + 1
    const index_type* column_ids,      // Size row_offsets[num_rows]
    const offset_type* row_bin_ids,
    offset_type actual_num_bins,
    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >& col_split_pattern)
{
    bool success =
        row_offsets &&
        column_ids &&
        row_bin_ids;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_to_col: Unacceptable input argument(s).");

        return false;
    }

    sparse_vectors<index_type, offset_type, offset_type> col_bins;

    sparse_vectors<const index_type, const offset_type, const offset_type> row_bins(
        num_rows, num_cols,
        row_offsets, column_ids,
        row_bin_ids);

    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >
        tmp_col_split_pattern;

    success =
        row_bins.get_transpose(
            col_bins)
        &&
        tmp_col_split_pattern.allocate(
            actual_num_bins, offset_type(1))
        &&
        split_pattern_to_bins(
            num_cols, num_rows,
            col_bins.vec_offsets(),
            col_bins.vec_ids(),
            col_bins.vec_values(),
            tmp_col_split_pattern.num_vecs(),
            tmp_col_split_pattern.vec_values());

    if(success)
    {
        col_split_pattern.swap(tmp_col_split_pattern);
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_binning_to_col: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_SPECTRAL_BINNING_H
