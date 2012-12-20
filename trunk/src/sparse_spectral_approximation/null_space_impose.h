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


#ifndef NULL_SPACE_IMPOSE_H
#define NULL_SPACE_IMPOSE_H

// -----------------------------------------------------------------------------

#include "sparse_algorithms/sparse_matrix_mult.h"
#include "sparse_vectors/sparse_vectors.h"
#include "math/vector_utils.h"
#include "math/precision_traits.h"
#include "math/complex_types.h"
#include "internal_api_error/internal_api_error.h"
#include <vector>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------

// sparsity structure of projected storage should already exist.

template
<
    typename index_type,
    typename value_type,
    typename index_type_2,
    typename offset_type_2
>
bool null_space_impose_project_residual(
    bool left,  // left or right
    index_type num_rows,
    index_type num_cols,
    index_type nullity,
    const value_type* basis,  // (left ? num_rows : num_cols) x nullity
    index_type basis_LD,
    const value_type* resid,  // (left ? num_cols : num_rows) x nullity
    index_type resid_LD,
    sparse_vectors<index_type_2, offset_type_2, value_type>& projected) // row based vectors
{
    bool success =
        num_rows == projected.num_vecs() &&
        num_cols == projected.max_size() &&
        nullity <= (left ? num_rows : num_cols) &&
        (left ? num_rows : num_cols) <= basis_LD &&
        (left ? num_cols : num_rows) <= resid_LD &&
        basis &&
        resid;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "null_space_impose_project_residual: Unacceptable input argument(s).");

        return false;
    }

    for(index_type i = 0; i < num_rows; ++i)
    {
        value_type* row_vals                = projected.vec_values_begin(i);
        const index_type_2* col_ids         = projected.vec_ids_begin(i);
        const index_type_2  num_vec_entries = projected.num_vec_entries(i);

        assert(row_vals);
        assert(col_ids);
        assert(num_vec_entries <= num_cols);

        for(index_type j = 0; j < num_vec_entries; ++j)
        {
            const index_type col = col_ids[j];

            value_type tmp = 0;
            const value_type* basis_it = basis + (left ? i : col);
            const value_type* resid_it = resid + (left ? col : i);

            for(index_type i_null = 0; i_null < nullity; ++i_null)
            {
                tmp += std::conj(*basis_it) * (*resid_it);

                basis_it += basis_LD;
                resid_it += resid_LD;
            }

            if(left)
                row_vals[j] = std::conj(tmp);
            else
                row_vals[j] = tmp;
        }
    }

    return success;
}

// -----------------------------------------------------------------------------

// FIXME: One could improve this method to be more efficient when working with
// Hermitian/skew-Hermitian/complex-symmetric matrices because the left and
// right null-spaces are related there.  However, this is not the most time
// consuming part of our main objective, so we let it be general and somewhat
// slower.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool null_space_impose(
// Input:
    index_type left_nullity,
    const value_type* left_basis,   // A.num_vecs() x left_nullity
    index_type left_basis_LD,
    index_type right_nullity,
    const value_type* right_basis,  // A.max_size() x right_nullity
    index_type right_basis_LD,

// Input/Output:
    sparse_vectors<index_type, offset_type, value_type>& A) // row based
{
    // A will be overwritten with the matrix closest to A in Frobenius
    // norm that also has the given left and right null-spaces.
    // We assume that basis for left and right null spaces are either
    // orthonormal or empty.

    const index_type num_rows = A.num_vecs();
    const index_type num_cols = A.max_size();
    const offset_type n_entries = A.num_entries();

    typedef typename precision_traits<value_type>::scalar scalar_type;

    const value_type* A_values = A.vec_values();

    scalar_type A_frob_norm_square = 0;

    for(offset_type i = 0; i < n_entries; ++i)
        A_frob_norm_square += std::abs_square(A_values[i]);

    const std::size_t max_iters = 1000; // MAGIC CONSTANT
    const scalar_type fuzz = 1;         // MAGIC CONSTANT

    const scalar_type sqtol = fuzz *
        A_frob_norm_square *
        std::numeric_limits<scalar_type>::epsilon() *
        std::numeric_limits<scalar_type>::epsilon() *
        scalar_type(std::max(num_rows, num_cols)) *
        scalar_type(std::max(num_rows, num_cols));

    std::vector<value_type> left_projected_store, right_projected_store;

    try
    {
        left_projected_store.resize(n_entries);
        right_projected_store.resize(n_entries);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("null_space_impose: Exception. ") + exc.what()));

        return false;
    }

    sparse_vectors<index_type, offset_type, value_type>
        left_projected(
            num_rows,
            num_cols,
            A.vec_offsets(),
            A.vec_ids(),
            &left_projected_store.front()),
        right_projected(
            num_rows,
            num_cols,
            A.vec_offsets(),
            A.vec_ids(),
            &right_projected_store.front());

    dense_vectors<index_type, value_type>
        left_lambda,  right_lambda,
        left_resid_1, right_resid_1,
        left_resid_2, right_resid_2;

    // Allocate local space for a specialized Uzawa CG iteration for this problem.

    const bool compute_lag_mult = false; // Lagrange multipliers.

    bool success =
        (compute_lag_mult ? left_lambda.allocate(left_nullity, num_cols) : true) &&
        left_resid_1.allocate (left_nullity, num_cols) &&
        left_resid_2.allocate (left_nullity, num_cols) &&
        (compute_lag_mult ? right_lambda.allocate(right_nullity, num_rows) : true) &&
        right_resid_1.allocate(right_nullity, num_rows) &&
        right_resid_2.allocate(right_nullity, num_rows);

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "null_space_impose_project_residual: Error 1.");

        return false;
    }

    if(success)
    {
        const value_type zero      = value_type(0);
        const value_type one       = value_type(1);
        const value_type minus_one = value_type(-1);

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4127 ) // conditional expression is constant
#endif
        if(compute_lag_mult)
#ifdef _MSC_VER
#pragma warning( pop )
#endif
        {
            left_lambda.fill(zero);
            right_lambda.fill(zero);
        }

        const dense_vectors<index_type, const value_type> left_basis_dv(
            left_nullity,
            num_rows,
            left_basis_LD,
            left_basis);

        const dense_vectors<index_type, const value_type> right_basis_dv(
            right_nullity,
            num_cols,
            right_basis_LD,
            right_basis);

        success =
            sparse_matrix_mult_trans(A, left_basis_dv, left_resid_1,  one) &&
            sparse_matrix_mult(A, right_basis_dv, right_resid_1, one) &&
            right_resid_2.axpby(right_resid_1, minus_one, zero) &&
            left_resid_2.axpby(left_resid_1, minus_one, zero);

        assert(success);

        std::size_t i_iter = 0;

        while(success && i_iter < max_iters)
        {
            const scalar_type R_norm_sq = right_resid_2.frobenius_norm_squared();
            const scalar_type L_norm_sq = left_resid_2.frobenius_norm_squared();

            if(R_norm_sq <= sqtol && L_norm_sq <= sqtol)
            {
                success = true;
                break;
            }

            success =
                null_space_impose_project_residual(
                    true, num_rows, num_cols,
                    left_nullity, left_basis, left_basis_LD, left_resid_1.vec_values(),
                    left_resid_1.leading_dimension(), left_projected)
                &&
                null_space_impose_project_residual(
                    false, num_rows, num_cols,
                    right_nullity, right_basis, right_basis_LD, right_resid_1.vec_values(),
                    right_resid_1.leading_dimension(), right_projected);

            assert(success);

            if(success)
            {
                // right one will now have the total
                vector_utils_add(
                    n_entries,
                    left_projected.vec_values(),
                    right_projected.vec_values(),
                    offset_type(1),
                    offset_type(1));

                const sparse_vectors<index_type, offset_type, value_type>&
                    total_projected = right_projected;

                const scalar_type proj_norm_sq = total_projected.frobenius_norm_squared();

                const value_type alpha = (R_norm_sq + L_norm_sq) / proj_norm_sq;

                vector_utils_axpby(
                    n_entries,
                    total_projected.vec_values(),
                    A.vec_values(),
                    -alpha,
                    value_type(1),
                    offset_type(1),
                    offset_type(1));

                success =
                    (compute_lag_mult ? left_lambda.axpby(left_resid_1, alpha, value_type(1)) : true) &&
                    (compute_lag_mult ? right_lambda.axpby(right_resid_1, alpha, value_type(1)) : true) &&
                    sparse_matrix_mult_trans(A, left_basis_dv,  left_resid_2, minus_one) &&
                    sparse_matrix_mult(A, right_basis_dv, right_resid_2, minus_one);

                assert(success);

                if(success)
                {
                    const scalar_type L_2_norm_sq = left_resid_2.frobenius_norm_squared();
                    const scalar_type R_2_norm_sq = right_resid_2.frobenius_norm_squared();

                    const value_type beta = (R_2_norm_sq + L_2_norm_sq)/(R_norm_sq + L_norm_sq);

                    success =
                        left_resid_1.axpby(left_resid_2, minus_one, beta) &&
                        right_resid_1.axpby(right_resid_2, minus_one, beta);

                    assert(success);
                }
            }

            ++i_iter;
        }
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "null_space_impose: Error 2.");
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
bool null_space_impose(
    const dense_vectors<index_type, value_type>& left_null_space,
    const dense_vectors<index_type, value_type>& right_null_space,
    sparse_vectors<index_type, offset_type, value_type>& A) // row based
{
    bool success =
        left_null_space.vec_size() == A.num_vecs() &&
        right_null_space.vec_size() == A.max_size() &&
        null_space_impose(
            left_null_space.num_vecs(),
            left_null_space.vec_values(),
            left_null_space.leading_dimension(),
            right_null_space.num_vecs(),
            right_null_space.vec_values(),
            right_null_space.leading_dimension(),
            A);

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "null_space_impose: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // NULL_SPACE_IMPOSE_H
