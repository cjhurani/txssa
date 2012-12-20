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


#ifndef DENSE_MATRIX_MULT_H
#define DENSE_MATRIX_MULT_H

// -----------------------------------------------------------------------------

#include "blas/blas_cpp_functions.h"
#include "blas/blas_char_check.h"
#include "platform/integral_type_range.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Output matrix will always be column-oriented.  One can play with transposes
// to achieve multiplication of both row and col-oriented inputs.
// -----------------------------------------------------------------------------

// This function exists solely to wrap BLAS gemm and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_mult(
                char  trans_a,
                char  trans_b,
          index_type  num_rows_C,
          index_type  num_cols_C,
          index_type  inner_size,
          value_type  alpha,
    const value_type* A_col_values,
          index_type  A_col_leading_dim,
    const value_type* B_col_values,
          index_type  B_col_leading_dim,
          value_type  beta,
          value_type* C_col_values,
          index_type  C_col_leading_dim)
{
    bool success =
        BLAS_char_check_trans(trans_a) &&
        BLAS_char_check_trans(trans_b) &&
        A_col_values &&
        B_col_values &&
        C_col_values &&
        num_rows_C <= C_col_leading_dim &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_rows_C) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_cols_C) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(inner_size) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(B_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(C_col_leading_dim);

    assert(success);

    if(success)
    {
        BLAS_int BLAS_m   = BLAS_int(num_rows_C);
        BLAS_int BLAS_n   = BLAS_int(num_cols_C);
        BLAS_int BLAS_k   = BLAS_int(inner_size);

        BLAS_int BLAS_lda = BLAS_int(A_col_leading_dim);
        BLAS_int BLAS_ldb = BLAS_int(B_col_leading_dim);
        BLAS_int BLAS_ldc = BLAS_int(C_col_leading_dim);

        BLAS_gemm(
            &trans_a, &trans_b,
            &BLAS_m, &BLAS_n, &BLAS_k,
            &alpha,
            const_cast<value_type*>(A_col_values), &BLAS_lda,
            const_cast<value_type*>(B_col_values), &BLAS_ldb,
            &beta,
            C_col_values, &BLAS_ldc);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_mult: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap BLAS syrk/herk and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_mult_herk(
                char  uplo,
                char  trans,
          index_type  matrix_size_C,
          index_type  inner_size,
          value_type  alpha,
    const value_type* A_col_values,
          index_type  A_col_leading_dim,
          value_type  beta,
          value_type* C_col_values,
          index_type  C_col_leading_dim)
{
    bool success =
        BLAS_char_check_uplo(uplo) &&
        BLAS_char_check_trans(trans) &&
        A_col_values &&
        C_col_values &&
        matrix_size_C <= C_col_leading_dim &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(matrix_size_C) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(inner_size) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(C_col_leading_dim);

    assert(success);

    if(success)
    {
        BLAS_int BLAS_n   = BLAS_int(matrix_size_C);
        BLAS_int BLAS_k   = BLAS_int(inner_size);

        BLAS_int BLAS_lda = BLAS_int(A_col_leading_dim);
        BLAS_int BLAS_ldc = BLAS_int(C_col_leading_dim);

        BLAS_herk(
            &uplo, &trans,
            &BLAS_n, &BLAS_k,
            &alpha,
            const_cast<value_type*>(A_col_values), &BLAS_lda,
            &beta,
            C_col_values, &BLAS_ldc);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_mult_herk: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap BLAS gemv and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_mult_vec(
                char  trans,
          index_type  num_rows,
          index_type  num_cols,
          value_type  alpha,
    const value_type* A_col_values,
          index_type  A_col_leading_dim,
    const value_type* X_col_values,
          index_type  X_inc,
          value_type  beta,
          value_type* Y_col_values,
          index_type  Y_inc)
{
    bool success =
        BLAS_char_check_trans(trans) &&
        A_col_values &&
        X_col_values &&
        Y_col_values &&
        0 < X_inc &&
        0 < Y_inc &&
        num_rows <= A_col_leading_dim &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(X_inc) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(Y_inc);

    assert(success);

    if(success)
    {
        BLAS_int BLAS_m    = BLAS_int(num_rows);
        BLAS_int BLAS_n    = BLAS_int(num_cols);

        BLAS_int BLAS_lda  = BLAS_int(A_col_leading_dim);
        BLAS_int BLAS_xinc = BLAS_int(X_inc);
        BLAS_int BLAS_yinc = BLAS_int(Y_inc);

        BLAS_gemv(
            &trans,
            &BLAS_m, &BLAS_n,
            &alpha,
            const_cast<value_type*>(A_col_values), &BLAS_lda,
            const_cast<value_type*>(X_col_values), &BLAS_xinc,
            &beta,
            Y_col_values, &BLAS_yinc);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_mult_vec: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap BLAS symm/hemm and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_mult_hemm(
                char  side,
                char  uplo,
          index_type  C_num_rows,
          index_type  C_num_cols,
          value_type  alpha,
    const value_type* A_col_values,
          index_type  A_col_leading_dim,
    const value_type* B_col_values,
          index_type  B_col_leading_dim,
          value_type  beta,
          value_type* C_col_values,
          index_type  C_col_leading_dim)
{
    bool success =
        BLAS_char_check_side(side) &&
        BLAS_char_check_uplo(uplo) &&
        A_col_values &&
        B_col_values &&
        C_col_values &&
        C_num_rows <= C_col_leading_dim &&
        C_num_rows <= B_col_leading_dim &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(C_num_rows) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(C_num_cols) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(B_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(C_col_leading_dim);

    assert(success);

    if(success)
    {
        BLAS_int BLAS_m   = BLAS_int(C_num_rows);
        BLAS_int BLAS_n   = BLAS_int(C_num_cols);

        BLAS_int BLAS_lda = BLAS_int(A_col_leading_dim);
        BLAS_int BLAS_ldb = BLAS_int(B_col_leading_dim);
        BLAS_int BLAS_ldc = BLAS_int(C_col_leading_dim);

        BLAS_hemm(
            &side, &uplo,
            &BLAS_m, &BLAS_n,
            &alpha,
            const_cast<value_type*>(A_col_values), &BLAS_lda,
            const_cast<value_type*>(B_col_values), &BLAS_ldb,
            &beta,
            C_col_values, &BLAS_ldc);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_mult_hemm: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_MULT_H
