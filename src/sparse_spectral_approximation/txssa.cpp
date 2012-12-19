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


// -----------------------------------------------------------------------------

#include "txssa.h"
#include "sparse_spectral_approximation/ssa_matrix_type.h"
#include "sparse_spectral_approximation/ssa_matrix_type_pinv_transpose.h"
#include "sparse_spectral_approximation/sparse_spectral_minimization.h"
#include "sparse_spectral_approximation/sparse_spectral_misfit_lhs_matrices.h"
#include "sparse_spectral_approximation/sparse_spectral_binning.h"
#include "sparse_vectors/sparse_vectors.h"
#include "p_norm_sparsity_matrix/p_norm_sparsity_dense_matrix.h"
#include "dense_algorithms/dense_matrix_utils.h"
#include "dense_vectors/dense_vectors.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <algorithm> // std::copy
#include <vector>
#include <stdexcept>
#include <complex>
#include <cassert>
#include <new>

// -----------------------------------------------------------------------------

#ifdef _MSC_VER
#pragma warning( disable : 4514 ) // unreferenced inline function has been removed
#pragma warning( disable : 4710 ) // function not inlined
#endif

// -----------------------------------------------------------------------------

const unsigned int ssa_impl_version = 0;

// -----------------------------------------------------------------------------

namespace
{
// -----------------------------------------------------------------------------

template<typename T>
void delete_catch(T* ptr, const char* func_name)
{
    try
    {
        delete ptr;
    }
    catch(...)
    {
        assert(false);
        internal_api_error_set_last(
            std::string(func_name) + ": Exception in deleting a pointer.");
    }
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool ssa_matrix_type_compute_AAT_from_ATA(
    index_type  size,
    const value_type* ATA_col_values,
    index_type  ATA_col_leading_dim,
    value_type* AAT_col_values,
    index_type  AAT_col_leading_dim,
    ssa_matrix_type type)
{
    bool success = false;

    // This is the only one implemented/needed now.
    if(type == ssa_matrix_type_complex_symmetric)
    {
        success = dense_matrix_utils_complex_sym_compute_AAT_from_ATA(
            size,
            ATA_col_values, ATA_col_leading_dim,
            AAT_col_values, AAT_col_leading_dim);
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_matrix_type_compute_AAT_from_ATA: Error");
    }

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
void ssa_B2TB2_chooser(
    const dense_vectors<index_type, value_type>& B1TB1,
    const dense_vectors<index_type, value_type>& B2TB2,
    ssa_matrix_type matrix_type,
    bool is_binned,
    const value_type*& tmp_B2TB2_col_values,
    typename precision_traits<value_type>::scalar& mult_factor)
{
    // The reason we compute B1B1T and then compute B2TB2 from it instead of
    // the other way around (this is when such the two matrices are related),
    // is that on some Intel MKL versions X'*X is slightly faster than X*X',
    // say, 5-10%.

    // hermitian implies normal, so only 6 combinations.

    // hermitian binned normal b2       mult_factor
    // 0         0      0      pass b2  2
    // 0         0      1      pass b1  2
    // 0         1      0      pass b2  2
    // 0         1      1      pass b1  2
    // 1         0      1      pass b1  2
    // 1         1      1      pass 0   1

    tmp_B2TB2_col_values = 0;

    const int is_hermitian = ssa_matrix_type_is_hermitian(matrix_type);
    const int is_normal    = ssa_matrix_type_is_normal(matrix_type);

    if(is_normal)
    {
        if(!is_hermitian || !is_binned)
        {
            tmp_B2TB2_col_values = B1TB1.vec_values();
        }
    }
    else
    {
        tmp_B2TB2_col_values = B2TB2.vec_values();
    }

    typedef typename precision_traits<value_type>::scalar scalar_type;

    mult_factor = tmp_B2TB2_col_values ? scalar_type(2) : scalar_type(1);
}

// -----------------------------------------------------------------------------

// User-given parameters for computing L_p norm based pattern.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool ssa_lpn_internal(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    typename precision_traits<value_type>::scalar sparsity_ratio,
    typename precision_traits<value_type>::scalar sparsity_norm_p,
    offset_type          max_num_bins,
    bool                 impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix)
{
    ssa_error_clear();

    const int is_abs_sym = ssa_matrix_type_is_abs_sym(matrix_type);

    bool success =
        ssa_matrix_type_undefined < matrix_type &&
        matrix_type < ssa_matrix_type_num_types &&
        (num_rows == num_cols || !is_abs_sym) &&
        col_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_lpn_internal: Unacceptable input argument(s).");

        return false;
    }

    dense_vectors<index_type, value_type> pinv_AT, left_null_space, right_null_space;

    dense_vectors<index_type, value_type>* left_null_space_ptr = 0;
    dense_vectors<index_type, value_type>* right_null_space_ptr = 0;

    if(impose_null_spaces)
    {
        left_null_space_ptr = &left_null_space;
        right_null_space_ptr = &right_null_space;
    }

    success =
        pinv_AT.allocate(
            num_cols, num_rows)
        &&
        dense_vectors_utils_copy(
            num_cols, num_rows,
            col_values, col_leading_dim,
            pinv_AT.vec_values(), pinv_AT.leading_dimension())
        &&
        ssa_matrix_type_pinv_transpose(
            num_rows, num_cols,
            pinv_AT.vec_values(), pinv_AT.leading_dimension(),
            matrix_type,
            left_null_space_ptr, right_null_space_ptr);

    if(success)
    {
        // ISSUE: Ideally, we should be using an equivalent of the near_zero_row_col
        // function written in the MATLAB code to compute min_num_nnz_per_row
        // and min_num_nnz_per_col.  We skip that step.

        const index_type min_num_nnz_per_row = right_null_space.num_vecs();
        const index_type min_num_nnz_per_col = left_null_space.num_vecs();

        std::vector< std::vector<index_type> > row_oriented_sparse_pat;

        if(is_abs_sym)
        {
            const int left_right_nullity_equal = ssa_matrix_type_is_left_right_nullity_equal(matrix_type);

            if(left_right_nullity_equal && min_num_nnz_per_row != min_num_nnz_per_col)
            {
                assert(false);

                internal_api_error_set_last(
                    "ssa_lpn_internal: Left and right nullity should be equal"
                    " because of matrix type, but not computed to be equal.");

                return false;
            }

            success = p_norm_sparsity_dense_matrix_abs_sym(
                sparsity_ratio, sparsity_norm_p,
                min_num_nnz_per_row,
                num_rows,
                col_values, col_leading_dim,
                row_oriented_sparse_pat);
        }
        else
        {
            success = p_norm_sparsity_dense_matrix_col_oriented(
                sparsity_ratio, sparsity_norm_p,
                min_num_nnz_per_row, min_num_nnz_per_col,
                num_rows, num_cols,
                col_values, col_leading_dim,
                row_oriented_sparse_pat);
        }

        if(success)
        {
            std::vector<index_type> size_per_row(num_rows);

            for(index_type row = 0; row < num_rows; ++row)
                size_per_row[row] = index_type(row_oriented_sparse_pat[row].size());

            sparse_vectors<index_type, offset_type, value_type>* out_mat_ptr =
                new (std::nothrow) sparse_vectors<index_type, offset_type, value_type>();

            if(!out_mat_ptr)
            {
                assert(false);

                internal_api_error_set_last(
                    "ssa_lpn_internal: Error in allocating matrix.");

                return false;
            }

            success = out_mat_ptr->allocate(num_rows, num_cols, &size_per_row.front());

            if(success)
            {
                for(index_type row = 0; row < num_rows; ++row)
                    std::copy(
                        row_oriented_sparse_pat[row].begin(),
                        row_oriented_sparse_pat[row].end(),
                        out_mat_ptr->vec_ids_begin(row));

                success =
                    ssa_internal(
                        num_rows, num_cols,
                        col_values, col_leading_dim,
                        out_mat_ptr->vec_offsets(),
                        out_mat_ptr->vec_ids(),
                        max_num_bins,
                        impose_null_spaces,
                        pinv_AT, left_null_space, right_null_space,
                        matrix_type,
                        out_mat_ptr->vec_values());

                if(success)
                {
                    out_matrix.row_offsets = out_mat_ptr->vec_offsets();
                    out_matrix.column_ids  = out_mat_ptr->vec_ids();
                    out_matrix.values      = out_mat_ptr->vec_values();
                    out_matrix.reserved    = out_mat_ptr;
                }
            }
        }
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_lpn_internal: Error");
    }

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool ssa_ids_internal(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    typename precision_traits<value_type>::scalar sparsity_ratio,
    typename precision_traits<value_type>::scalar sparsity_norm_p,
    index_type           min_num_nnz_per_row,
    index_type           min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix)
{
    ssa_error_clear();

    const int is_abs_sym = ssa_matrix_type_is_abs_sym(matrix_type);

    bool success =
        ssa_matrix_type_undefined < matrix_type &&
        matrix_type < ssa_matrix_type_num_types &&
        (num_rows == num_cols || !is_abs_sym) &&
        (min_num_nnz_per_row == min_num_nnz_per_col || !is_abs_sym) &&
        col_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_ids_internal: Unacceptable input argument(s).");

        return false;
    }

    std::vector< std::vector<index_type> > row_oriented_sparse_pat;

    if(is_abs_sym)
    {
        success = p_norm_sparsity_dense_matrix_abs_sym(
            sparsity_ratio, sparsity_norm_p,
            min_num_nnz_per_row,
            num_rows,
            col_values, col_leading_dim,
            row_oriented_sparse_pat);
    }
    else
    {
        success = p_norm_sparsity_dense_matrix_col_oriented(
            sparsity_ratio, sparsity_norm_p,
            min_num_nnz_per_row, min_num_nnz_per_col,
            num_rows, num_cols,
            col_values, col_leading_dim,
            row_oriented_sparse_pat);
    }

    if(success)
    {
        std::vector<index_type> size_per_row(num_rows);

        for(index_type row = 0; row < num_rows; ++row)
            size_per_row[row] = index_type(row_oriented_sparse_pat[row].size());

        sparse_vectors<index_type, offset_type, value_type>* out_mat_ptr =
            new (std::nothrow) sparse_vectors<index_type, offset_type, value_type>();

        if(!out_mat_ptr)
        {
            assert(false);

            internal_api_error_set_last(
                "ssa_ids_internal: Error in allocating matrix.");

            return false;
        }

        success = out_mat_ptr->allocate(num_rows, num_cols, &size_per_row.front());

        if(success)
        {
            for(index_type row = 0; row < num_rows; ++row)
                std::copy(
                    row_oriented_sparse_pat[row].begin(),
                    row_oriented_sparse_pat[row].end(),
                    out_mat_ptr->vec_ids_begin(row));

            out_matrix.row_offsets = out_mat_ptr->vec_offsets();
            out_matrix.column_ids  = out_mat_ptr->vec_ids();
            out_matrix.values      = out_mat_ptr->vec_values();
            out_matrix.reserved    = out_mat_ptr;
        }
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_ids_internal: Error");
    }

    return success;
}

// -----------------------------------------------------------------------------

// For real matrices
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool ssa_internal(
    index_type num_rows,
    index_type num_cols,
    const value_type* col_values,      // num_rows x num_cols
    index_type col_leading_dim,
    const offset_type* row_offsets,    // num_rows + 1
    const index_type* column_ids,      // row_offsets[num_rows]
    offset_type max_num_bins,
    bool impose_null_spaces,
    const dense_vectors<index_type, value_type>& pinv_AT,
    const dense_vectors<index_type, value_type>& left_null_space,
    const dense_vectors<index_type, value_type>& right_null_space,
    ssa_matrix_type matrix_type,
    value_type* out_row_values)  // row_offsets[num_rows]
{
    const int is_square = (num_rows == num_cols);
    const int is_hermitian = ssa_matrix_type_is_hermitian(matrix_type);
    const int is_normal    = ssa_matrix_type_is_normal(matrix_type);

    bool success =
        !(
            is_hermitian ||
            is_normal
        ) ||
        is_square;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_internal: Unacceptable input argument(s) in real version.");

        return false;
    }

    std::vector<offset_type> row_bin_ids;

    try
    {
        row_bin_ids.resize(row_offsets[num_rows]);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("ssa_internal: Exception in the real version. ") + exc.what()));

        return false;
    }

    dense_vectors<index_type, value_type> B1TB1, B2TB2;

    offset_type actual_num_bins;

    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >
        row_split_pattern, col_split_pattern;

    success =
        sparse_spectral_binning_row(
            num_rows, num_cols,
            col_values, col_leading_dim,
            row_offsets, column_ids,
            max_num_bins, actual_num_bins,
            row_split_pattern,
            row_bin_ids.size() ? &row_bin_ids.front() : 0)
        &&
        // If hermitian, don't have to compute col_split_pattern
        (is_hermitian ? true : sparse_spectral_binning_to_col(
            num_rows, num_cols,
            row_offsets, column_ids,
            row_bin_ids.size() ? &row_bin_ids.front() : 0,
            actual_num_bins,
            col_split_pattern))
        &&
        B1TB1.allocate(
            num_cols, num_cols)
        &&
        // If normal, don't have to allocate or compute B2TB2
        (is_normal ? true : B2TB2.allocate(
            num_rows, num_rows))
        &&
        sparse_spectral_misfit_lhs_matrices(
            num_rows,
            num_cols,
            pinv_AT.vec_values(), pinv_AT.leading_dimension(),
            B2TB2.vec_values(),   B2TB2.leading_dimension(),
            B1TB1.vec_values(),   B1TB1.leading_dimension());

    if(success)
    {
        const value_type* tmp_B2TB2_col_values = 0;
        value_type mult_factor;

        const bool is_binned = max_num_bins != 0;

        ssa_B2TB2_chooser<index_type, value_type>(
            B1TB1,
            B2TB2,
            matrix_type,
            is_binned,
            tmp_B2TB2_col_values,
            mult_factor);

        success = sparse_spectral_minimization(
            num_rows, num_cols,
            row_offsets, column_ids,
            actual_num_bins,
            row_bin_ids.size() ? &row_bin_ids.front() : 0,
            row_split_pattern.vec_values(),
            // If hermitian, reuse row_split_pattern as col_split_pattern
            is_hermitian ? row_split_pattern.vec_values() : col_split_pattern.vec_values(),
            impose_null_spaces,
            tmp_B2TB2_col_values, num_rows,
            B1TB1.vec_values(), B1TB1.leading_dimension(),
            pinv_AT.vec_values(), pinv_AT.leading_dimension(),
            left_null_space,
            right_null_space,
            out_row_values,
            mult_factor);
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_internal: Error in real version.");
    }

    return success;
}

// -----------------------------------------------------------------------------

// For complex matrices
template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
bool ssa_internal(
    index_type num_rows,
    index_type num_cols,
    const std::complex<scalar_type>* col_values,    // num_rows x num_cols
    index_type col_leading_dim,
    const offset_type* row_offsets,    // num_rows + 1
    const index_type* column_ids,      // row_offsets[num_rows]
    offset_type max_num_bins,
    bool impose_null_spaces,
    const dense_vectors<index_type, std::complex<scalar_type> >& pinv_AT,
    const dense_vectors<index_type, std::complex<scalar_type> >& left_null_space,
    const dense_vectors<index_type, std::complex<scalar_type> >& right_null_space,
    ssa_matrix_type matrix_type,
    std::complex<scalar_type>* out_row_values)  // row_offsets[num_rows]
{
    const int is_square = (num_rows == num_cols);
    const int is_hermitian               = ssa_matrix_type_is_hermitian(matrix_type);
    const int is_normal                  = ssa_matrix_type_is_normal(matrix_type);
    const int is_real_part_symmetric     = ssa_matrix_type_is_real_part_symmetric(matrix_type);
    const int is_imag_part_symmetric     = ssa_matrix_type_is_imag_part_symmetric(matrix_type);
    const int is_AAT_computable_from_ATA = ssa_matrix_type_is_AAT_computable_from_ATA(matrix_type);

    bool success =
        !(
            is_hermitian ||
            is_normal ||
            is_real_part_symmetric ||
            is_imag_part_symmetric ||
            is_AAT_computable_from_ATA
        ) ||
        is_square;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_internal: Unacceptable input argument(s) in complex version.");

        return false;
    }

    std::vector<offset_type> real_row_bin_ids, imag_row_bin_ids;

    try
    {
        real_row_bin_ids.resize(row_offsets[num_rows]);
        imag_row_bin_ids.resize(row_offsets[num_rows]);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("ssa_internal: Exception in complex version. ") + exc.what()));

        return false;
    }

    dense_vectors<index_type, std::complex<scalar_type> > B1TB1, B2TB2;

    offset_type real_actual_num_bins, imag_actual_num_bins;

    dense_vectors<offset_type, sparse_vectors_ids<index_type, offset_type> >
        real_row_split_pattern, real_col_split_pattern,
        imag_row_split_pattern, imag_col_split_pattern;

    success =
        sparse_spectral_binning_row(
            num_rows, num_cols,
            col_values, col_leading_dim,
            row_offsets, column_ids,
            max_num_bins,
            real_actual_num_bins, imag_actual_num_bins,
            real_row_split_pattern, imag_row_split_pattern,
            real_row_bin_ids.size() ? &real_row_bin_ids.front() : 0,
            imag_row_bin_ids.size() ? &imag_row_bin_ids.front() : 0)
        &&
        // If real part hermitian, don't have to compute real_col_split_pattern
        (is_real_part_symmetric ? true : sparse_spectral_binning_to_col(
            num_rows, num_cols,
            row_offsets, column_ids,
            real_row_bin_ids.size() ? &real_row_bin_ids.front() : 0,
            real_actual_num_bins,
            real_col_split_pattern))
        &&
        // If imag part hermitian, don't have to compute imag_col_split_pattern
        (is_imag_part_symmetric ? true : sparse_spectral_binning_to_col(
            num_rows, num_cols,
            row_offsets, column_ids,
            imag_row_bin_ids.size() ? &imag_row_bin_ids.front() : 0,
            imag_actual_num_bins,
            imag_col_split_pattern))
        &&
        B1TB1.allocate(
            num_cols, num_cols)
        &&
        // If normal, don't have to allocate or compute B2TB2
        (is_normal ? true : B2TB2.allocate(
            num_rows, num_rows));

    if(success)
    {
        if(is_AAT_computable_from_ATA)
        {
            success =
                sparse_spectral_misfit_lhs_matrices(
                    num_rows,
                    num_cols,
                    pinv_AT.vec_values(), pinv_AT.leading_dimension(),
                    reinterpret_cast<std::complex<scalar_type>*>(0), num_rows,
                    B1TB1.vec_values(),   B1TB1.leading_dimension())
                &&
                ssa_matrix_type_compute_AAT_from_ATA(
                    num_rows,
                    B1TB1.vec_values(), B1TB1.leading_dimension(),
                    B2TB2.vec_values(), B2TB2.leading_dimension(),
                    matrix_type);
        }
        else
        {
            success = sparse_spectral_misfit_lhs_matrices(
                num_rows,
                num_cols,
                pinv_AT.vec_values(), pinv_AT.leading_dimension(),
                B2TB2.vec_values(),   B2TB2.leading_dimension(),
                B1TB1.vec_values(),   B1TB1.leading_dimension());
        }

        if(success)
        {
            const std::complex<scalar_type>* tmp_B2TB2_col_values = 0;
            scalar_type mult_factor;

            const bool is_binned = max_num_bins != 0;

            ssa_B2TB2_chooser<index_type, std::complex<scalar_type> >(
                B1TB1,
                B2TB2,
                matrix_type,
                is_binned,
                tmp_B2TB2_col_values,
                mult_factor);

            success = sparse_spectral_minimization(
                num_rows, num_cols,
                row_offsets, column_ids,
                real_actual_num_bins, imag_actual_num_bins,
                real_row_bin_ids.size() ? &real_row_bin_ids.front() : 0,
                imag_row_bin_ids.size() ? &imag_row_bin_ids.front() : 0,
                real_row_split_pattern.vec_values(), imag_row_split_pattern.vec_values(),
                is_real_part_symmetric ? real_row_split_pattern.vec_values() : real_col_split_pattern.vec_values(),
                is_imag_part_symmetric ? imag_row_split_pattern.vec_values() : imag_col_split_pattern.vec_values(),
                impose_null_spaces,
                tmp_B2TB2_col_values, num_rows,
                B1TB1.vec_values(), B1TB1.leading_dimension(),
                pinv_AT.vec_values(), pinv_AT.leading_dimension(),
                left_null_space,
                right_null_space,
                out_row_values,
                mult_factor);
        }
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_internal: Error in complex version.");
    }

    return success;
}

// -----------------------------------------------------------------------------

} // namespace

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
ssa_csr<index_type, offset_type, value_type>::ssa_csr()
    :
        row_offsets(0),
        column_ids(0),
        values(0),
        reserved(0)
{
}

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
ssa_csr<index_type, offset_type, value_type>::~ssa_csr()
{
    delete_catch(reinterpret_cast
        <const sparse_vectors<index_type, offset_type, value_type>*>(
            reserved), "ssa_csr::~ssa_csr()");

    reserved = 0;
}

// -----------------------------------------------------------------------------

// User-given pattern.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
int ssa_pat(
    index_type            num_rows,
    index_type            num_cols,
    const value_type*     col_values,
    index_type            col_leading_dim,
    const offset_type*    row_offsets,
    const index_type*     column_ids,
    offset_type           max_num_bins,
    bool                  impose_null_spaces,
    enum ssa_matrix_type  matrix_type,
    value_type*           out_row_values)
{
    ssa_error_clear();

    bool success =
        ssa_matrix_type_undefined < matrix_type &&
        matrix_type < ssa_matrix_type_num_types &&
        (num_rows == num_cols || !ssa_matrix_type_is_abs_sym(matrix_type)) &&
        col_values &&
        row_offsets &&
        column_ids &&
        out_row_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_pat: Unacceptable input argument(s).");

        return false;
    }

    dense_vectors<index_type, value_type>
        pinv_AT, left_null_space, right_null_space;

    dense_vectors<index_type, value_type>* left_null_space_ptr = 0;
    dense_vectors<index_type, value_type>* right_null_space_ptr = 0;

    if(impose_null_spaces)
    {
        left_null_space_ptr = &left_null_space;
        right_null_space_ptr = &right_null_space;
    }

    success =
        pinv_AT.allocate(
            num_cols, num_rows)
        &&
        dense_vectors_utils_copy(
            num_cols, num_rows,
            col_values, col_leading_dim,
            pinv_AT.vec_values(), pinv_AT.leading_dimension())
        &&
        ssa_matrix_type_pinv_transpose(
            num_rows, num_cols,
            pinv_AT.vec_values(), pinv_AT.leading_dimension(),
            matrix_type,
            left_null_space_ptr, right_null_space_ptr)
        &&
        ssa_internal(
            num_rows, num_cols,
            col_values, col_leading_dim,
            row_offsets, column_ids,
            max_num_bins,
            impose_null_spaces,
            pinv_AT, left_null_space, right_null_space,
            matrix_type,
            out_row_values);

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_pat: Error.");
    }

    return success ? 0 : 1;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
int ssa_lpn(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    value_type           sparsity_ratio,
    value_type           sparsity_norm_p,
    offset_type          max_num_bins,
    bool                 impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix)
{
    bool success = ssa_lpn_internal<index_type, offset_type, value_type>(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces,
        matrix_type,
        out_matrix);

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_lpn: Error in real version");
    }

    return success ? 0 : 1;
}

template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
int ssa_lpn(
    index_type                       num_rows,
    index_type                       num_cols,
    const std::complex<scalar_type>* col_values,
    index_type                       col_leading_dim,
    scalar_type                      sparsity_ratio,
    scalar_type                      sparsity_norm_p,
    offset_type                      max_num_bins,
    bool                             impose_null_spaces,
    enum ssa_matrix_type             matrix_type,
    ssa_csr<index_type, offset_type, std::complex<scalar_type> >& out_matrix)
{
    bool success = ssa_lpn_internal<index_type, offset_type, std::complex<scalar_type> >(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces,
        matrix_type,
        out_matrix);

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_lpn: Error in complex version");
    }

    return success ? 0 : 1;
}


// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
int ssa_ids(
    index_type           num_rows,
    index_type           num_cols,
    const value_type*    col_values,
    index_type           col_leading_dim,
    value_type           sparsity_ratio,
    value_type           sparsity_norm_p,
    index_type           min_num_nnz_per_row,
    index_type           min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, value_type>& out_matrix)
{
    bool success = ssa_ids_internal<index_type, offset_type, value_type>(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        out_matrix);

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_ids: Error in real version");
    }

    return success ? 0 : 1;
}

template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
int ssa_ids(
    index_type           num_rows,
    index_type           num_cols,
    const std::complex<scalar_type>* col_values,
    index_type           col_leading_dim,
    scalar_type          sparsity_ratio,
    scalar_type          sparsity_norm_p,
    index_type           min_num_nnz_per_row,
    index_type           min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    ssa_csr<index_type, offset_type, std::complex<scalar_type> >& out_matrix)
{
    bool success = ssa_ids_internal<index_type, offset_type, std::complex<scalar_type> >(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        out_matrix);

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "ssa_ids: Error in complex version");
    }

    return success ? 0 : 1;
}

// -----------------------------------------------------------------------------

int ssa_error_size()
{
    return internal_api_error_size();
}

int ssa_error_string(int i, const char** ptr_to_error_string)
{
    return internal_api_error_string(i, ptr_to_error_string);
}

int ssa_error_clear()
{
    return internal_api_error_clear();
}

// -----------------------------------------------------------------------------

/* User-given pattern */
int ssa_d_pat(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    double*              out_row_values)
{
    int ret = ssa_pat(
        num_rows, num_cols,
        col_values, col_leading_dim,
        row_offsets, column_ids,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        out_row_values);

    return ret;
}

/* User-given parameters for computing L_p norm based pattern. */
int ssa_d_lpn(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_d_csr*    out_matrix)
{
    ssa_csr<int, int, double>* csr = new (std::nothrow) ssa_csr<int, int, double>;

    int ret = ssa_lpn(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = csr->values;
    out_matrix->reserved    = csr;

    return ret;
}

int ssa_d_ids(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_d_csr*    out_matrix)
{
    ssa_csr<int, int, double>* csr = new (std::nothrow) ssa_csr<int, int, double>;

    int ret = ssa_ids(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = csr->values;
    out_matrix->reserved    = csr;

    return ret;
}

// -----------------------------------------------------------------------------

int ssa_s_pat(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    float*               out_row_values)
{
    int ret = ssa_pat(
        num_rows, num_cols,
        col_values, col_leading_dim,
        row_offsets, column_ids,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        out_row_values);

    return ret;
}

int ssa_s_lpn(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_s_csr*    out_matrix)
{
    ssa_csr<int, int, float>* csr = new (std::nothrow) ssa_csr<int, int, float>;

    int ret = ssa_lpn(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = csr->values;
    out_matrix->reserved    = csr;

    return ret;
}

int ssa_s_ids(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_s_csr*    out_matrix)
{
    ssa_csr<int, int, float>* csr = new (std::nothrow) ssa_csr<int, int, float>;

    int ret = ssa_ids(
        num_rows, num_cols,
        col_values, col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = csr->values;
    out_matrix->reserved    = csr;

    return ret;
}

// -----------------------------------------------------------------------------

int ssa_z_pat(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    double*              out_row_values)
{
    int ret = ssa_pat<int,int,std::complex<double> >(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<double>*>(col_values), col_leading_dim,
        row_offsets, column_ids,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        reinterpret_cast<std::complex<double>*>(out_row_values));

    return ret;
}

int ssa_z_lpn(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_z_csr*    out_matrix)
{
    ssa_csr<int, int, std::complex<double> >* csr = new (std::nothrow) ssa_csr<int, int, std::complex<double> >;

    int ret = ssa_lpn(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<double>*>(col_values), col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = reinterpret_cast<double*>(csr->values);
    out_matrix->reserved    = csr;

    return ret;
}

int ssa_z_ids(
    int                  num_rows,
    int                  num_cols,
    const double*        col_values,
    int                  col_leading_dim,
    double               sparsity_ratio,
    double               sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_z_csr*    out_matrix)
{
    ssa_csr<int, int, std::complex<double> >* csr = new (std::nothrow) ssa_csr<int, int, std::complex<double> >;

    int ret = ssa_ids(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<double>*>(col_values), col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = reinterpret_cast<double*>(csr->values);
    out_matrix->reserved    = csr;

    return ret;
}

// -----------------------------------------------------------------------------

int ssa_c_pat(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    const int*           row_offsets,
    const int*           column_ids,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    float*               out_row_values)
{
    int ret = ssa_pat<int,int,std::complex<float> >(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<float>*>(col_values), col_leading_dim,
        row_offsets, column_ids,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        reinterpret_cast<std::complex<float>*>(out_row_values));

    return ret;
}

int ssa_c_lpn(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  max_num_bins,
    int                  impose_null_spaces,
    enum ssa_matrix_type matrix_type,
    struct ssa_c_csr*    out_matrix)
{
    ssa_csr<int, int, std::complex<float> >* csr = new (std::nothrow) ssa_csr<int, int, std::complex<float> >;

    int ret = ssa_lpn(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<float>*>(col_values), col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        max_num_bins,
        impose_null_spaces == 0 ? false : true,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = reinterpret_cast<float*>(csr->values);
    out_matrix->reserved    = csr;

    return ret;
}

int ssa_c_ids(
    int                  num_rows,
    int                  num_cols,
    const float*         col_values,
    int                  col_leading_dim,
    float                sparsity_ratio,
    float                sparsity_norm_p,
    int                  min_num_nnz_per_row,
    int                  min_num_nnz_per_col,
    enum ssa_matrix_type matrix_type,
    struct ssa_c_csr*    out_matrix)
{
    ssa_csr<int, int, std::complex<float> >* csr = new (std::nothrow) ssa_csr<int, int, std::complex<float> >;

    int ret = ssa_ids(
        num_rows, num_cols,
        reinterpret_cast<const std::complex<float>*>(col_values), col_leading_dim,
        sparsity_ratio, sparsity_norm_p,
        min_num_nnz_per_row, min_num_nnz_per_col,
        matrix_type,
        *csr);

    out_matrix->row_offsets = csr->row_offsets;
    out_matrix->column_ids  = csr->column_ids;
    out_matrix->values      = reinterpret_cast<float*>(csr->values);
    out_matrix->reserved    = csr;

    return ret;
}

// -----------------------------------------------------------------------------

void ssa_d_csr_deallocate(struct ssa_d_csr* matrix)
{
    if(matrix)
    {
        matrix->row_offsets = 0;
        matrix->column_ids = 0;
        matrix->values = 0;

        delete_catch(
            reinterpret_cast<
                const ssa_csr<int, int, double>*>(
                    matrix->reserved), "ssa_d_csr_deallocate");

        matrix->reserved = 0;
    }
}

void ssa_s_csr_deallocate(struct ssa_s_csr* matrix)
{
    if(matrix)
    {
        matrix->row_offsets = 0;
        matrix->column_ids = 0;
        matrix->values = 0;

        delete_catch(
            reinterpret_cast<
                const ssa_csr<int, int, float>*>(
                    matrix->reserved), "ssa_s_csr_deallocate");

        matrix->reserved = 0;
    }
}

void ssa_z_csr_deallocate(struct ssa_z_csr* matrix)
{
    if(matrix)
    {
        matrix->row_offsets = 0;
        matrix->column_ids = 0;
        matrix->values = 0;

        delete_catch(
            reinterpret_cast<
                const ssa_csr<int, int, std::complex<double> >*>(
                    matrix->reserved), "ssa_z_csr_deallocate");

        matrix->reserved = 0;
    }
}

void ssa_c_csr_deallocate(struct ssa_c_csr* matrix)
{
    if(matrix)
    {
        matrix->row_offsets = 0;
        matrix->column_ids = 0;
        matrix->values = 0;

        delete_catch(
            reinterpret_cast<
                const ssa_csr<int, int, std::complex<float> >*>(
                    matrix->reserved), "ssa_c_csr_deallocate");

        matrix->reserved = 0;
    }
}

// -----------------------------------------------------------------------------
// instantiate the template functions explicitly.

#define SSA_INSTANTIATE_PAT_LPN(index, offset, scalar)         \
                                                               \
template TXSSA_API int ssa_pat<index, offset, scalar>(         \
    index            num_rows,                                 \
    index            num_cols,                                 \
    const scalar*    col_values,                               \
    index            col_leading_dim,                          \
    const offset*    row_offsets,                              \
    const index*     column_ids,                               \
    offset           max_num_bins,                             \
    bool             impose_null_spaces,                       \
    enum ssa_matrix_type matrix_type,                          \
    scalar*          out_row_values);                          \
                                                               \
template TXSSA_API int ssa_lpn<index, offset, scalar>(         \
    index           num_rows,                                  \
    index           num_cols,                                  \
    const scalar*   col_values,                                \
    index           col_leading_dim,                           \
    scalar          sparsity_ratio,                            \
    scalar          sparsity_norm_p,                           \
    offset          max_num_bins,                              \
    bool            impose_null_spaces,                        \
    enum ssa_matrix_type matrix_type,                          \
    ssa_csr<index, offset, scalar>& out_matrix);               \
                                                               \
template TXSSA_API int ssa_pat<index, offset, std::complex<scalar> >(    \
    index            num_rows,                                 \
    index            num_cols,                                 \
    const std::complex<scalar>* col_values,                    \
    index            col_leading_dim,                          \
    const offset*    row_offsets,                              \
    const index*     column_ids,                               \
    offset           max_num_bins,                             \
    bool             impose_null_spaces,                       \
    enum ssa_matrix_type matrix_type,                          \
    std::complex<scalar>* out_row_values);                     \
                                                               \
template TXSSA_API int ssa_lpn<index, offset, scalar>(         \
    index           num_rows,                                  \
    index           num_cols,                                  \
    const std::complex<scalar>* col_values,                    \
    index           col_leading_dim,                           \
    scalar          sparsity_ratio,                            \
    scalar          sparsity_norm_p,                           \
    offset          max_num_bins,                              \
    bool            impose_null_spaces,                        \
    enum ssa_matrix_type matrix_type,                          \
    ssa_csr<index, offset, std::complex<scalar> >& out_matrix)

#define SSA_INSTANTIATE_IDS(index, offset, scalar)             \
template TXSSA_API int ssa_ids<index, offset, scalar>(         \
    index           num_rows,                                  \
    index           num_cols,                                  \
    const scalar*   col_values,                                \
    index           col_leading_dim,                           \
    scalar          sparsity_ratio,                            \
    scalar          sparsity_norm_p,                           \
    index           min_num_nnz_per_row,                       \
    index           min_num_nnz_per_col,                       \
    enum ssa_matrix_type matrix_type,                          \
    ssa_csr<index, offset, scalar>& out_matrix);               \
template TXSSA_API int ssa_ids<index, offset, scalar>(         \
    index           num_rows,                                  \
    index           num_cols,                                  \
    const std::complex<scalar>* col_values,                    \
    index           col_leading_dim,                           \
    scalar          sparsity_ratio,                            \
    scalar          sparsity_norm_p,                           \
    index           min_num_nnz_per_row,                       \
    index           min_num_nnz_per_col,                       \
    enum ssa_matrix_type matrix_type,                          \
    ssa_csr<index, offset, std::complex<scalar> >& out_matrix)

#define SSA_INSTANTIATE_CSR(index, offset, scalar)                       \
    template ssa_csr<index, offset, scalar>::ssa_csr();                  \
    template ssa_csr<index, offset, scalar>::~ssa_csr();                 \
    template ssa_csr<index, offset, std::complex<scalar> >::ssa_csr();   \
    template ssa_csr<index, offset, std::complex<scalar> >::~ssa_csr()

#define SSA_INSTANTIATE_PAT_LPN_IDS_CSR(index, offset, scalar)            \
        SSA_INSTANTIATE_PAT_LPN(index, offset, scalar);                   \
        SSA_INSTANTIATE_PAT_LPN(unsigned index, unsigned offset, scalar); \
        SSA_INSTANTIATE_IDS(index, offset, scalar);                       \
        SSA_INSTANTIATE_IDS(unsigned index, unsigned offset, scalar);     \
        SSA_INSTANTIATE_CSR(index, offset, scalar);                       \
        SSA_INSTANTIATE_CSR(unsigned index, unsigned offset, scalar)

#define SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(index, offset)          \
        SSA_INSTANTIATE_PAT_LPN_IDS_CSR(index, offset, float);            \
        SSA_INSTANTIATE_PAT_LPN_IDS_CSR(index, offset, double)

// Commenting out the less used ones.

//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(char, char);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(char, short);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(char, int);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(char, long);


//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(short, short);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(short, int);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(short, long);

SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(int, int);
//SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(int, long);

#if defined(__x86_64__) 
SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(long, long);
#endif

#if defined(_WIN64)
SSA_INSTANTIATE_PAT_LPN_IDS_CSR_all_float(long long, long long);
#endif

