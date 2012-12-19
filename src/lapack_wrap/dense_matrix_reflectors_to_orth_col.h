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


#ifndef DENSE_MATRIX_REFLECTORS_TO_ORTH_COL_H
#define DENSE_MATRIX_REFLECTORS_TO_ORTH_COL_H

// -----------------------------------------------------------------------------

#include "lapack/lapack_cpp_functions.h"
#include "platform/integral_type_range.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <cmath>      // std::real
#include <cstddef>
#include <limits>
#include <cassert>

// -----------------------------------------------------------------------------

// Call this function to get lwork for ORGQR/UNGQR.

template<typename index_type, typename value_type>
std::size_t dense_matrix_reflectors_to_orth_col_lwork(
    index_type num_rows,
    index_type num_cols,
    index_type num_reflectors,
    index_type A_col_leading_dim,
    value_type /**/)
{
    bool success =
        num_rows <= A_col_leading_dim &&
        num_cols <= num_rows &&
        num_reflectors <= num_cols &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_reflectors) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim);

    assert(success);

    value_type work = value_type();

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_ref = LAPACK_int(num_reflectors);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(-1);

        int info = 0;

        LAPACK_ungqr(
            &LAPACK_m, &LAPACK_n, &LAPACK_ref,
            0, &LAPACK_lda,
            0,
            &work, &LAPACK_wsz,
            &info);

        typedef typename precision_traits<value_type>::scalar scalar_type;
        success =
            (scalar_type(std::size_t(std::real(work))) == std::real(work)) &&
            info == 0;

        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_reflectors_to_orth_col_lwork: Error.");

    return success ? std::size_t(std::real(work)) : std::numeric_limits<std::size_t>::max();
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap LAPACK OR/UNGQR and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_reflectors_to_orth_col(
    index_type  num_rows,
    index_type  num_cols,
    index_type  num_reflectors,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    const value_type* tau_reflectors,
    value_type* work,
    std::size_t work_size)
{
    bool success =
        A_col_values &&
        (num_reflectors == 0 || tau_reflectors) &&
        work &&
        num_rows <= A_col_leading_dim &&
        num_cols <= num_rows &&
        num_reflectors <= num_cols &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_reflectors) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(work_size);

    assert(success);

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_ref = LAPACK_int(num_reflectors);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(work_size);

        int info = 0;

        LAPACK_ungqr(
            &LAPACK_m, &LAPACK_n, &LAPACK_ref,
            A_col_values, &LAPACK_lda,
            const_cast<value_type*>(tau_reflectors),
            work, &LAPACK_wsz,
            &info);

        success = (info == 0);
        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_reflectors_to_orth_col: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_REFLECTORS_TO_ORTH_COL_H
