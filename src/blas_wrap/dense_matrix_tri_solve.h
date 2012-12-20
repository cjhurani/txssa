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


#ifndef DENSE_MATRIX_TRI_SOLVE_H
#define DENSE_MATRIX_TRI_SOLVE_H

// -----------------------------------------------------------------------------

#include "blas/blas_cpp_functions.h"
#include "blas/blas_char_check.h"
#include "platform/integral_type_range.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------

// This function exists solely to wrap BLAS trsm and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_tri_solve(
                char  side,
                char  uplo,
                char  trans_a,
                char  diag,
          index_type  num_rows_B,
          index_type  num_cols_B,
          value_type  alpha,
    const value_type* A_col_values,
          index_type  A_col_leading_dim,
          value_type* B_col_values,
          index_type  B_col_leading_dim)
{
    bool success =
        BLAS_char_check_trans(trans_a) &&
        BLAS_char_check_side(side) &&
        BLAS_char_check_uplo(uplo) &&
        BLAS_char_check_diag(diag) &&
        A_col_values &&
        B_col_values &&
        num_rows_B <= B_col_leading_dim &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_rows_B) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(num_cols_B) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<BLAS_int>::in_non_negative_range(B_col_leading_dim);

    assert(success);

    if(success)
    {
        BLAS_int BLAS_m   = BLAS_int(num_rows_B);
        BLAS_int BLAS_n   = BLAS_int(num_cols_B);
        BLAS_int BLAS_lda = BLAS_int(A_col_leading_dim);
        BLAS_int BLAS_ldb = BLAS_int(B_col_leading_dim);

        BLAS_trsm(
            &side, &uplo, &trans_a, &diag,
            &BLAS_m, &BLAS_n,
            &alpha,
            const_cast<value_type*>(A_col_values), &BLAS_lda,
            B_col_values, &BLAS_ldb);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_tri_solve: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_TRI_SOLVE_H
