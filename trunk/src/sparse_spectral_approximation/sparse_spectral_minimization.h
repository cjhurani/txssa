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


#ifndef SPARSE_SPECTRAL_MINIMIZATION_H
#define SPARSE_SPECTRAL_MINIMIZATION_H

// -----------------------------------------------------------------------------

#include "sparse_spectral_approximation/sparse_spectral_misfit_lhs.h"
#include "sparse_spectral_approximation/sparse_spectral_misfit_rhs.h"
#include "sparse_spectral_approximation/null_space_impose.h"
#include "sparse_vectors/sparse_vectors.h"
#include "lapack_wrap/dense_matrix_linear_hpd.h"
#include "dense_vectors/dense_vectors.h"
#include "math/vector_utils.h"
#include "internal_api_error/internal_api_error.h"
#include <complex>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute the sparse spectral approximation given the input matrix
// and data defining the minimization problem.
// -----------------------------------------------------------------------------

// mult_factor is used to multiply the solution obtained.  It is useful in that
// it allows us to not have to multiply RHS_col_values, which is a matrix and
// the resulting product will consume extra space.

// For real matrices.
template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_spectral_minimization(
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets,    // Size num_rows + 1
    const index_type* column_ids,      // Size row_offsets[num_rows]
    offset_type actual_num_bins,
    const offset_type* row_bin_values,   // Size row_offsets[num_rows]
    const sparse_vectors_ids<index_type, offset_type>* row_split_pattern, // actual_num_bins
    const sparse_vectors_ids<index_type, offset_type>* col_split_pattern, // actual_num_bins
    bool impose_null_spaces,
    const value_type* B2TB2_col_values,  // num_rows x num_rows
    index_type  B2TB2_col_leading_dim,
    const value_type* B1B1T_col_values,  // num_cols x num_cols
    index_type  B1B1T_col_leading_dim,
    const value_type* RHS_col_values,    // num_rows x num_cols
    index_type RHS_col_leading_dim,
    const dense_vectors<index_type, value_type>& left_null_space,
    const dense_vectors<index_type, value_type>& right_null_space,
    value_type* out_row_values,  // Size row_offsets[num_rows]
    value_type mult_factor)
{
    bool success = out_row_values && row_bin_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_minimization: Unacceptable input argument(s).");

        return false;
    }

    dense_vectors<offset_type, value_type> LS_A, LS_b; // LS = Least squares

    success =
        LS_A.allocate(
            actual_num_bins, actual_num_bins)
        &&
        LS_b.allocate(
            index_type(1), actual_num_bins)
        &&
        sparse_spectral_misfit_lhs(
            num_rows, num_cols,
            B2TB2_col_values, B2TB2_col_leading_dim,
            B1B1T_col_values, B1B1T_col_leading_dim,
            actual_num_bins,
            row_split_pattern,
            col_split_pattern,
            LS_A.vec_values(),
            LS_A.leading_dimension())
        &&
        sparse_spectral_misfit_rhs(
            num_rows, num_cols,
            RHS_col_values, RHS_col_leading_dim,
            actual_num_bins,
            col_split_pattern,
            LS_b.vec_values())
        &&
        dense_matrix_linear_hpd_solve(
            'U',
            actual_num_bins,
            offset_type(1),
            LS_A.vec_values(),
            LS_A.leading_dimension(),
            LS_b.vec_values(),
            actual_num_bins);

    if(success)
    {
        value_type* LS_b_vals = LS_b.vec_values();

        if(mult_factor != 1)
            vector_utils_axpby(
                LS_b.max_size(),
                reinterpret_cast<const value_type*>(0),
                LS_b_vals,
                value_type(0),
                mult_factor,
                offset_type(1),
                offset_type(1));

        for(offset_type i = 0; i < row_offsets[num_rows]; ++i)
        {
            out_row_values[i] = LS_b_vals[row_bin_values[i]];
        }

        if((left_null_space.num_vecs() || right_null_space.num_vecs()) && impose_null_spaces)
        {
            sparse_vectors<index_type, offset_type, value_type> approximation(
                num_rows, num_cols, const_cast<offset_type*>(row_offsets), const_cast<index_type*>(column_ids),
                out_row_values);

            success = null_space_impose(
                left_null_space,
                right_null_space,
                approximation);
        }
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_minimization: Error.");
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
bool sparse_spectral_minimization(
    index_type num_rows,
    index_type num_cols,
    const offset_type* row_offsets,    // Size num_rows + 1
    const index_type* column_ids,      // Size row_offsets[num_rows]
    offset_type real_actual_num_bins,
    offset_type imag_actual_num_bins,
    const offset_type* real_row_bin_values,   // Size row_offsets[num_rows]
    const offset_type* imag_row_bin_values,   // Size row_offsets[num_rows]
    const sparse_vectors_ids<index_type, offset_type>* real_row_split_pattern, // real_actual_num_bins
    const sparse_vectors_ids<index_type, offset_type>* imag_row_split_pattern, // imag_actual_num_bins
    const sparse_vectors_ids<index_type, offset_type>* real_col_split_pattern, // real_actual_num_bins
    const sparse_vectors_ids<index_type, offset_type>* imag_col_split_pattern, // imag_actual_num_bins
    bool impose_null_spaces,
    const std::complex<scalar_type>* B2TB2_col_values,  // num_rows x num_rows
    index_type  B2TB2_col_leading_dim,
    const std::complex<scalar_type>* B1B1T_col_values,  // num_cols x num_cols
    index_type  B1B1T_col_leading_dim,
    const std::complex<scalar_type>* RHS_col_values,    // num_rows x num_cols
    index_type RHS_col_leading_dim,
    const dense_vectors<index_type, std::complex<scalar_type> >& left_null_space,
    const dense_vectors<index_type, std::complex<scalar_type> >& right_null_space,
    std::complex<scalar_type>* out_row_values,  // Size row_offsets[num_rows]
    scalar_type mult_factor)
{
    bool success =
        out_row_values &&
        real_row_bin_values &&
        imag_row_bin_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_minimization: Unacceptable input argument(s).");

        return false;
    }

    dense_vectors<offset_type, scalar_type> LS_A, LS_b; // LS = Least squares

    const offset_type actual_num_bins =
        real_actual_num_bins + imag_actual_num_bins;

    success =
        LS_A.allocate(
            actual_num_bins, actual_num_bins)
        &&
        LS_b.allocate(
            index_type(1), actual_num_bins)
        &&
        sparse_spectral_misfit_lhs(
            num_rows, num_cols,
            B2TB2_col_values, B2TB2_col_leading_dim,
            B1B1T_col_values, B1B1T_col_leading_dim,
            real_actual_num_bins, imag_actual_num_bins,
            real_row_split_pattern,
            imag_row_split_pattern,
            real_col_split_pattern,
            imag_col_split_pattern,
            LS_A.vec_values(),
            LS_A.leading_dimension())
        &&
        sparse_spectral_misfit_rhs(
            num_rows, num_cols,
            RHS_col_values, RHS_col_leading_dim,
            real_actual_num_bins, imag_actual_num_bins,
            real_col_split_pattern,
            imag_col_split_pattern,
            LS_b.vec_values())
        &&
        dense_matrix_linear_hpd_solve(
            'U',
            actual_num_bins,
            offset_type(1),
            LS_A.vec_values(),
            LS_A.leading_dimension(),
            LS_b.vec_values(),
            LS_b.leading_dimension());

    if(success)
    {
        scalar_type* LS_b_vals = LS_b.vec_values();

        if(mult_factor != 1)
            vector_utils_axpby(
                LS_b.max_size(),
                reinterpret_cast<const scalar_type*>(0),
                LS_b_vals,
                scalar_type(0),
                mult_factor,
                offset_type(1),
                offset_type(1));

        typedef std::complex<scalar_type> value_type;

        for(offset_type i = 0; i < row_offsets[num_rows]; ++i)
        {
            out_row_values[i] = value_type(
                LS_b_vals[real_row_bin_values[i]],
                LS_b_vals[imag_row_bin_values[i] + real_actual_num_bins]);
        }

        if((left_null_space.num_vecs() || right_null_space.num_vecs()) && impose_null_spaces)
        {
            sparse_vectors<index_type, offset_type, value_type> approximation(
                num_rows, num_cols, const_cast<offset_type*>(row_offsets), const_cast<index_type*>(column_ids),
                out_row_values);

            success = null_space_impose(
                left_null_space,
                right_null_space,
                approximation);
        }
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_minimization: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_SPECTRAL_MINIMIZATION_H
