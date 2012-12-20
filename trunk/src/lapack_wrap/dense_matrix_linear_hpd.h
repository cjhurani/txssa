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


#ifndef DENSE_MATRIX_LINEAR_HPD_H
#define DENSE_MATRIX_LINEAR_HPD_H

// -----------------------------------------------------------------------------

#include "blas_wrap/dense_matrix_tri_solve.h"
#include "lapack/lapack_cpp_functions.h"
#include "blas/blas_char_check.h"
#include "platform/integral_type_range.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Function related to solving linear SPD systems.
// -----------------------------------------------------------------------------

// This function exists solely to wrap LAPACK POSV and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_linear_hpd_solve(
    char uplo,
    index_type  matrix_size,
    index_type  num_rhs,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type* B_col_values,
    index_type  B_col_leading_dim)
{
    bool success =
        BLAS_char_check_uplo(uplo) &&
        A_col_values &&
        B_col_values &&
        matrix_size <= A_col_leading_dim &&
        matrix_size <= B_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(matrix_size) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rhs) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(B_col_leading_dim);

    assert(success);

    if(success)
    {
        LAPACK_int LAPACK_n    = LAPACK_int(matrix_size);
        LAPACK_int LAPACK_nrhs = LAPACK_int(num_rhs);
        LAPACK_int LAPACK_lda  = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_ldb  = LAPACK_int(B_col_leading_dim);

        int info = 0;
        LAPACK_posv(
            &uplo, &LAPACK_n, &LAPACK_nrhs,
            A_col_values, &LAPACK_lda,
            B_col_values, &LAPACK_ldb,
            &info);

        success = (info == 0);
        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_linear_hpd_solve: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap LAPACK POTRF and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_linear_hpd_factor(
    char uplo,
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim)
{
    bool success =
        BLAS_char_check_uplo(uplo) &&
        A_col_values &&
        matrix_size <= A_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(matrix_size) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim);

    assert(success);

    if(success)
    {
        LAPACK_int LAPACK_n   = LAPACK_int(matrix_size);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);

        int info = 0;
        LAPACK_potrf(
            &uplo, &LAPACK_n,
            A_col_values, &LAPACK_lda,
            &info);

        success = (info == 0);
        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_linear_hpd_factor: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Solve X * A = B for X, when A is SPD, and X, A, B are stored column-wise.
// X will be output in space for B.  A will be over-written.
// This function is needed because LAPACK doesn't provide it.
// See http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=1530 (for LU
// talk on "How to solve XA = B").
// Solving A * X' = B' doesn't work in POTRS.

template<typename index_type, typename value_type>
bool dense_matrix_linear_hpd_solve_flip(
    char uplo,
    index_type  A_size,
    index_type  X_num_rows,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type* B_col_values,
    index_type  B_col_leading_dim)
{
    bool success =
        BLAS_char_check_uplo(uplo) &&
        A_col_values &&
        B_col_values &&
        A_size <= A_col_leading_dim &&
        X_num_rows <= B_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_size) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(X_num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(B_col_leading_dim);

    assert(success);

    if(success)
    {
        const int uplo_is_up = BLAS_char_uplo_is_up(uplo);
        const char trans_1 = uplo_is_up ? 'N' : 'C';
        const char trans_2 = uplo_is_up ? 'C' : 'N';
        const char unit_diag = 'N';

        success =
            dense_matrix_linear_hpd_factor(
                uplo, A_size, A_col_values, A_col_leading_dim)
            &&
            dense_matrix_tri_solve(
                'R', uplo, trans_1, unit_diag,
                X_num_rows, A_size, value_type(1),
                A_col_values, A_col_leading_dim,
                B_col_values, B_col_leading_dim)
                &&
            dense_matrix_tri_solve(
                'R', uplo, trans_2, unit_diag,
                X_num_rows, A_size, value_type(1),
                A_col_values, A_col_leading_dim,
                B_col_values, B_col_leading_dim);

        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_linear_hpd_solve_flip: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_LINEAR_HPD_H
