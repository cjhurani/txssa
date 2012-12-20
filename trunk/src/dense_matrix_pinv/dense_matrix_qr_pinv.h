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


#ifndef DENSE_MATRIX_QR_PINV_H
#define DENSE_MATRIX_QR_PINV_H

// -----------------------------------------------------------------------------

#include "dense_algorithms/dense_matrix_permute.h"
#include "dense_algorithms/dense_matrix_utils.h"
#include "dense_vectors/dense_vectors.h"
#include "dense_vectors/dense_vectors_utils.h"
#include "lapack_wrap/dense_matrix_linear_hpd.h"
#include "lapack_wrap/dense_matrix_tri_invert.h"
#include "lapack_wrap/dense_matrix_QR.h"
#include "lapack_wrap/dense_matrix_reflectors_mult.h"
#include "blas_wrap/dense_matrix_mult.h"
#include "math/precision_traits.h"
#include "platform/integral_type_range.h"
#include "internal_api_error/internal_api_error.h"
#include <stdexcept>
#include <vector>
#include <cstddef>
#include <cassert>

// -----------------------------------------------------------------------------

// dense_matrix_qr_pinv_horzcat_U_fullrank_B_flat: Compute pinv of [U B]
//
// horzcat means the full matrix is flat rectangular concatenation of U and B.
// U_fullrank means the leftmost square is upper triangular and full rank.
// B_flat means the remainder is flat rectangular (more columns than rows).
// This function should be used when num_rows <= B_num_cols.  It will work
// correctly even if that's not the case.  But it will be less efficient
// than dense_matrix_qr_pinv_horzcat_U_fullrank_B_tall.

template<typename index_type, typename value_type>
bool dense_matrix_qr_pinv_horzcat_U_fullrank_B_flat(
    index_type  num_rows,
    index_type  B_num_cols,
    value_type* U_col_values,           // num_rows x num_rows
    index_type  U_col_leading_dim,
    value_type* B_col_values,           // num_rows x B_num_cols
    index_type  B_col_leading_dim,
    value_type* inv_U_B_col_values,     // num_rows x B_num_cols
    index_type  inv_U_B_col_leading_dim,
    value_type* work)                   // num_rows x num_rows at least if B_num_cols > 0
{
    // On output:
    // U_col_values will contain top num_rows x num_rows part of pinv( [U, B] )
    // B_col_values will contain transpose of bottom B_num_cols x num_rows part of pinv( [U, B] )
    // inv_U_B_col_values will contain U^-1 B.

    bool success =
        (  num_rows == 0 || U_col_values) &&
        (B_num_cols == 0 || work) &&
        (B_num_cols == 0 || B_col_values) &&
        (B_num_cols == 0 || inv_U_B_col_values) &&
        num_rows <= U_col_leading_dim &&
        num_rows <= B_col_leading_dim &&
        num_rows <= inv_U_B_col_leading_dim
        &&
        // Invert U in-place.
        dense_matrix_tri_invert(
            'U', 'N', num_rows, U_col_values, U_col_leading_dim)
        &&
        // Fill strict lower part of U with zeros.  We will use
        // full U space below.
        dense_matrix_utils_fill_strict_lower(
            num_rows, num_rows, U_col_values, U_col_leading_dim, value_type(0));

    assert(success);

    if(success && 0 < B_num_cols)
    {
        const char CC_uplo = 'U';

        success =
            // C = inv_U * B;
            dense_matrix_mult(
                'N', 'N',
                num_rows, B_num_cols, num_rows, value_type(1),
                U_col_values, U_col_leading_dim,
                B_col_values, B_col_leading_dim,
                value_type(0), inv_U_B_col_values, inv_U_B_col_leading_dim)
            &&
            // Form C * C' and place it in "work"
            dense_matrix_mult_herk(
                CC_uplo, 'N',
                num_rows, B_num_cols, value_type(1),
                inv_U_B_col_values, inv_U_B_col_leading_dim, value_type(0),
                work, num_rows)
            &&
            // Add identity to work.
            dense_matrix_utils_diagonal_add(
                num_rows, num_rows,
                work,
                num_rows, value_type(1))
            &&
            // Replace inv(U) with inv(work) * inv(U).
            dense_matrix_linear_hpd_solve(
                CC_uplo, num_rows, num_rows,
                work, num_rows, U_col_values, U_col_leading_dim)
            &&
            // Compute the matrix to be put in B_col_values.
            dense_matrix_mult(
                'C', 'N',
                num_rows, B_num_cols, num_rows,
                value_type(1), U_col_values, U_col_leading_dim,
                inv_U_B_col_values, inv_U_B_col_leading_dim,
                value_type(0), B_col_values, B_col_leading_dim);

        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_qr_pinv_horzcat_U_fullrank_B_flat: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// dense_matrix_qr_pinv_horzcat_U_fullrank_B_tall: Compute pinv of [U B]
//
// horzcat means the full matrix is flat rectangular concatenation of U and B.
// U_fullrank means the leftmost square is upper triangular and full rank.
// B_tall means the remainder is tall rectangular (more rows than columns).
// This function should be used when num_rows >= B_num_cols.  It will work
// correctly even if that's not the case.  But it will be less efficient
// than dense_matrix_qr_pinv_horzcat_U_fullrank_B_flat.

// This one uses Woodbury identity.

template<typename index_type, typename value_type>
bool dense_matrix_qr_pinv_horzcat_U_fullrank_B_tall(
    index_type  num_rows,
    index_type  B_num_cols,
    value_type* U_col_values,           // num_rows x num_rows
    index_type  U_col_leading_dim,
    value_type* B_col_values,           // num_rows x B_num_cols
    index_type  B_col_leading_dim,
    value_type* inv_U_B_col_values,     // num_rows x B_num_cols
    index_type  inv_U_B_col_leading_dim,
    value_type* work)                   // B_num_cols x B_num_cols at least
{
    // On output:
    // U_col_values will contain top num_rows x num_rows part of pinv( [U, B] )
    // B_col_values will contain transpose of bottom B_num_cols x num_rows part of pinv( [U, B] )
    // inv_U_B_col_values will contain U^-1 B.

    bool success =
        (  num_rows == 0 || U_col_values) &&
        (B_num_cols == 0 || work) &&
        (B_num_cols == 0 || B_col_values) &&
        (B_num_cols == 0 || inv_U_B_col_values) &&
        num_rows <= U_col_leading_dim &&
        num_rows <= B_col_leading_dim &&
        num_rows <= inv_U_B_col_leading_dim
        &&
        // Invert U in-place.
        dense_matrix_tri_invert(
            'U', 'N', num_rows, U_col_values, U_col_leading_dim)
        &&
        // Fill strict lower part of U with zeros.
        dense_matrix_utils_fill_strict_lower(
            num_rows, num_rows,
            U_col_values, U_col_leading_dim,
            value_type(0));

    assert(success);

    if(success && 0 < B_num_cols)
    {
        const char CC_uplo = 'U';

        success =
            // C = inv_U * B
            // inv_U_B_col_values contains C
            dense_matrix_mult(
                'N', 'N',
                num_rows, B_num_cols, num_rows, value_type(1),
                U_col_values, U_col_leading_dim,
                B_col_values, B_col_leading_dim,
                value_type(0), inv_U_B_col_values, inv_U_B_col_leading_dim)
            &&
            // "B space" = inv(U)' * C
            // B_col_values contains inv(U)' * C
            dense_matrix_mult(
                'C', 'N',
                num_rows, B_num_cols, num_rows, value_type(1),
                U_col_values, U_col_leading_dim,
                inv_U_B_col_values, inv_U_B_col_leading_dim,
                value_type(0), B_col_values, B_col_leading_dim)
            &&
            // Form C' * C and place it in "work"
            dense_matrix_mult_herk(
                CC_uplo, 'C',
                B_num_cols, num_rows, value_type(1),
                inv_U_B_col_values, inv_U_B_col_leading_dim, value_type(0),
                work, B_num_cols)
            &&
            // Add identity to work.
            dense_matrix_utils_diagonal_add(
                B_num_cols, B_num_cols,
                work, B_num_cols,
                value_type(1))
            &&
            dense_matrix_linear_hpd_solve_flip(
                CC_uplo, B_num_cols, num_rows,
                work, B_num_cols,
                B_col_values, B_col_leading_dim)
            &&
            // B_col_values now contains
            // inv(U)' * C * inv(C'*C + I)
            // inv(U) <- (I - C * inv(C'*C + I) * C') * inv(U)
            dense_matrix_mult(
                'N', 'C',
                num_rows, num_rows, B_num_cols,
                value_type(-1), inv_U_B_col_values, inv_U_B_col_leading_dim,
                B_col_values, B_col_leading_dim,
                value_type(1), U_col_values, U_col_leading_dim);

        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_qr_pinv_horzcat_U_fullrank_B_tall: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Automatically choose between tall and flat.

template<typename index_type, typename value_type>
bool dense_matrix_qr_pinv_horzcat_U_fullrank_B_auto(
    index_type  num_rows,
    index_type  B_num_cols,
    value_type* U_col_values,           // num_rows x num_rows
    index_type  U_col_leading_dim,
    value_type* B_col_values,           // num_rows x B_num_cols
    index_type  B_col_leading_dim,
    value_type* inv_U_B_col_values,     // num_rows x B_num_cols
    index_type  inv_U_B_col_leading_dim,
    value_type* work)                   // min(num_rows, B_num_cols) square.
{
    return num_rows <= B_num_cols ?
        dense_matrix_qr_pinv_horzcat_U_fullrank_B_flat(
            num_rows,
            B_num_cols,
            U_col_values,
            U_col_leading_dim,
            B_col_values,
            B_col_leading_dim,
            inv_U_B_col_values,
            inv_U_B_col_leading_dim,
            work)
            :
        dense_matrix_qr_pinv_horzcat_U_fullrank_B_tall(
            num_rows,
            B_num_cols,
            U_col_values,
            U_col_leading_dim,
            B_col_values,
            B_col_leading_dim,
            inv_U_B_col_values,
            inv_U_B_col_leading_dim,
            work);
}

// -----------------------------------------------------------------------------

// If rnull is not null, it will contain right null-space as output.
// Caller deallocates that.  A_col_values will be overwritten with pinv(A)'.

template<typename index_type, typename value_type>
bool dense_matrix_qr_pinv_transpose(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    dense_vectors<index_type, value_type>* lnull,
    dense_vectors<index_type, value_type>* rnull)
{
    bool success = false;

    if(num_rows == 0 || num_cols == 0)
    {
        success =
            rnull->use_memory(0, num_cols, index_type(1), 0) &&
            lnull->use_memory(0, num_rows, index_type(1), 0);

        if(!success)
            internal_api_error_set_last(
                "dense_matrix_qr_pinv_transpose: Error 1.");

        assert(success);

        return success;
    }

    success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    if(!success)
    {
        assert(false);

        if(!success)
            internal_api_error_set_last(
                "dense_matrix_qr_pinv_transpose: Unacceptable input argument(s).");

        return false;
    }

    typedef typename precision_traits<value_type>::scalar precision_scalar;

    const std::size_t lwork = dense_matrix_QR_pivoted_lwork<index_type, value_type>(
        num_rows, num_cols, A_col_leading_dim);

    const std::size_t rwork_size = dense_matrix_QR_pivoted_rwork_size(num_cols, value_type());

    const index_type min_rows_cols = std::min(num_rows, num_cols);

    std::vector<index_type> pivots;
    std::vector<value_type> tau_reflectors;
    std::vector<value_type> work;
    std::vector<precision_scalar> rwork;
    std::vector<value_type> Q_work;

    try
    {
        pivots.resize(num_cols);
        tau_reflectors.resize(min_rows_cols);
        work.resize(lwork);
        rwork.resize(rwork_size);
        Q_work.resize(std::size_t(num_rows) * std::size_t(num_cols));
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("dense_matrix_qr_pinv_transpose: Exception 1. ") + exc.what()));

        return false;
    }

    precision_scalar* rwork_ptr = rwork.size() > 0 ? &rwork.front() : 0;

    success =
        dense_matrix_QR_pivoted(
            num_rows, num_cols, A_col_values, A_col_leading_dim,
            &pivots.front(), &tau_reflectors.front(),
            &work.front(), lwork, rwork_ptr)
            &&
        dense_matrix_utils_copy_strict_lower(
            num_rows, num_cols,
            A_col_values, A_col_leading_dim,
            &Q_work.front(), num_rows);

    if(!success)
    {
        assert(false);

        if(!success)
            internal_api_error_set_last(
                "dense_matrix_qr_pinv_transpose: Error 2.");

        return false;
    }

    const precision_scalar fuzz = 100; // MAGIC CONSTANT

    const index_type right_null_size = dense_matrix_QR_pivoted_right_null_size(
        num_rows, num_cols, A_col_values, A_col_leading_dim, fuzz);

    assert(right_null_size <= num_cols);

    const index_type U_size = index_type(num_cols - right_null_size);
    const index_type B_num_cols = right_null_size;

    value_type* U_col_values = A_col_values;
    const index_type  U_col_leading_dim = A_col_leading_dim;

    value_type* B_col_values = A_col_values + 
        std::size_t(A_col_leading_dim) * std::size_t(U_size);
    const index_type B_col_leading_dim = A_col_leading_dim;

    const index_type rnull_num_rows = num_cols;

    dense_vectors<index_type, value_type> rnull_tmp, rnull_tmp_perm;

    success =
        rnull_tmp.allocate(right_null_size, rnull_num_rows) &&
        rnull_tmp_perm.allocate(right_null_size, rnull_num_rows);

    if(!success)
    {
        assert(false);

        if(!success)
            internal_api_error_set_last(
                "dense_matrix_qr_pinv_transpose: Error 3.");

        return false;
    }

    std::vector<value_type> rect_pinv_work;

    try
    {
        if(U_size <= B_num_cols)
        {
            rect_pinv_work.resize(std::size_t(num_rows) * std::size_t(num_rows));
        }
        else
        {
            rect_pinv_work.resize(std::size_t(B_num_cols) * std::size_t(B_num_cols));
        }
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("dense_matrix_qr_pinv_transpose: Exception 2. ") + exc.what()));

        return false;
    }

    success = dense_matrix_qr_pinv_horzcat_U_fullrank_B_auto(
        U_size, B_num_cols,
        U_col_values, U_col_leading_dim,
        B_col_values, B_col_leading_dim,
        rnull_tmp.vec_values(),
        rnull_tmp.leading_dimension(),
        B_num_cols ? &rect_pinv_work.front() : 0);

    if(success && U_size < num_cols && rnull)
    {
        // Fill lower rnull_tmp space

        const index_type diff = index_type(num_cols - U_size);

        success =
            dense_vectors_utils_fill(
                diff, diff,
                rnull_tmp.vec_values() + U_size,
                rnull_tmp.leading_dimension(), value_type(0))
            &&
            dense_matrix_utils_diagonal_add(
                diff, diff,
                rnull_tmp.vec_values() + U_size,
                rnull_tmp.leading_dimension(), value_type(-1));
    }

    if(success && rnull)
    {
        if(rnull_tmp.num_vecs())
        {
            success =
                dense_matrix_permute_rows(
                    rnull_tmp.vec_size(),
                    rnull_tmp.num_vecs(),
                    rnull_tmp.vec_values(),
                    rnull_tmp.leading_dimension(),
                    &pivots.front(),
                    index_type(1), // 1 since coming from LAPACK.
                    rnull_tmp_perm.vec_values(),
                    rnull_tmp_perm.leading_dimension())
                &&
                // Give an orthogonal basis for the null space.
                dense_matrix_QR_orth_col_space_for_full_rank(
                    rnull_tmp_perm.vec_size(),
                    rnull_tmp_perm.num_vecs(),
                    rnull_tmp_perm.vec_values(),
                    rnull_tmp_perm.leading_dimension());

            if(!success)
            {
                assert(false);

                if(!success)
                    internal_api_error_set_last(
                        "dense_matrix_qr_pinv_transpose: Error 4.");

                return false;
            }
        }

        rnull->swap(rnull_tmp_perm);
    }

    success =
        dense_matrix_utils_transpose_in_place(
            U_size, U_col_values, U_col_leading_dim)
        &&
        dense_vectors_utils_fill(
            num_cols,
            index_type(num_rows - U_size),
            A_col_values + U_size,
            A_col_leading_dim,
            value_type(0));

    if(!success)
    {
        assert(false);

        if(!success)
            internal_api_error_set_last(
                "dense_matrix_qr_pinv_transpose: Error 5.");

        return false;
    }

    const char mqr_side  = 'L';
    const char mqr_trans = 'N';

    const std::size_t mqr_work_size = dense_matrix_reflectors_mult_lwork<index_type, value_type>(
        mqr_side, mqr_trans,
        num_rows, num_cols, min_rows_cols,
        num_rows, A_col_leading_dim);

    std::vector<value_type> reflectors_mult_work;

    try
    {
        reflectors_mult_work.resize(mqr_work_size);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("dense_matrix_qr_pinv_transpose: Exception 3. ") + exc.what()));

        return false;
    }

    success = dense_matrix_reflectors_mult(
        mqr_side, mqr_trans,
        num_rows, num_cols, min_rows_cols,
        &Q_work.front(), num_rows,
        &tau_reflectors.front(),
        A_col_values, A_col_leading_dim,
        &reflectors_mult_work.front(), mqr_work_size);

    if(success)
    {
        if(lnull && right_null_size + num_rows >= num_cols)
        {
            dense_vectors<index_type, value_type> lnull_tmp;

            const index_type left_null_size =
                index_type(index_type(right_null_size + num_rows) - num_cols);

            success =
                lnull_tmp.allocate(
                    left_null_size, num_rows)
                &&
                lnull_tmp.fill(
                    value_type(0))
                &&
                dense_matrix_utils_diagonal_add(
                    num_rows, left_null_size,
                    lnull_tmp.vec_values() + index_type(num_rows - left_null_size),
                    lnull_tmp.leading_dimension(), value_type(1))
                &&
                dense_matrix_reflectors_mult(
                    mqr_side, mqr_trans,
                    num_rows, left_null_size, min_rows_cols,
                    &Q_work.front(), num_rows,
                    &tau_reflectors.front(),
                    lnull_tmp.vec_values(), lnull_tmp.leading_dimension(),
                    &reflectors_mult_work.front(), mqr_work_size);

            if(!success)
            {
                assert(false);

                if(!success)
                    internal_api_error_set_last(
                        "dense_matrix_qr_pinv_transpose: Error 6.");

                return false;
            }

            lnull->swap(lnull_tmp);
        }

        success =
            dense_matrix_permute_cols(
                num_rows, num_cols,
                A_col_values, A_col_leading_dim,
                &pivots.front(), index_type(1), // 1 since coming from LAPACK.
                &Q_work.front(), num_rows)
            &&
            dense_vectors_utils_copy(
                num_cols,
                num_rows,
                &Q_work.front(), num_rows,
                A_col_values, A_col_leading_dim);

        if(!success)
        {
            assert(false);

            if(!success)
                internal_api_error_set_last(
                    "dense_matrix_qr_pinv_transpose: Error 7.");

            return false;
        }
    }

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_qr_pinv_transpose: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_QR_PINV_H
