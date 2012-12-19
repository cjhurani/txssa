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


#ifndef DENSE_MATRIX_QR_H
#define DENSE_MATRIX_QR_H

// -----------------------------------------------------------------------------

#include "lapack_wrap/dense_matrix_reflectors_to_orth_col.h"
#include "lapack/lapack_cpp_functions.h"
#include "math/precision_traits.h"
#include "platform/integral_type_range.h"
#include "internal_api_error/internal_api_error.h"
#include "cpp/std_utils.h"
#include <cassert>
#include <cstddef>
#include <vector>
#include <limits>
#include <stdexcept>
#include <complex>
#include <cmath>      // std::real
#include <algorithm>  // std::{min, max}

// -----------------------------------------------------------------------------
// Functions related to pivoted and unpivoted QR factorization.
// -----------------------------------------------------------------------------

// Call one of these functions to get rwork size for GEQP3.

template<typename index_type, typename value_type>
std::size_t dense_matrix_QR_pivoted_rwork_size(index_type /*num_cols*/, value_type /**/)
{
    return 0;
}

template<typename index_type, typename value_type>
std::size_t dense_matrix_QR_pivoted_rwork_size(index_type num_cols, std::complex<value_type> /**/)
{
    return 2 * std::size_t(num_cols);
}

// -----------------------------------------------------------------------------

// Call this function to get lwork for GEQP3.

template<typename index_type, typename value_type>
std::size_t dense_matrix_QR_pivoted_lwork(
    index_type num_rows,
    index_type num_cols,
    index_type A_col_leading_dim)
{
    bool success =
        num_rows <= A_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim);

    assert(success);

    value_type work = value_type();

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(-1);

        int info = 0;
        LAPACK_geqp3(
            &LAPACK_m, &LAPACK_n,
            0, &LAPACK_lda,
            0, 0,
            &work, &LAPACK_wsz, 0,
            &info);

        typedef typename precision_traits<value_type>::scalar scalar_type;
        success = 
            (scalar_type(std::size_t(std::real(work))) == std::real(work)) &&
            info == 0;

        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_QR_pivoted_lwork: Error.");

    return success ? std::size_t(std::real(work)) : std::numeric_limits<std::size_t>::max();
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap LAPACK GEQP3 and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_QR_pivoted(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    index_type* pivots,
    value_type* tau_reflectors,
    value_type* work,
    std::size_t work_size,
    typename precision_traits<value_type>::scalar* rwork)
{
    // rwork can be null.

    bool success =
        A_col_values &&
        pivots &&
        (std::min(num_rows, num_cols) == 0 || tau_reflectors) &&
        work &&
        num_rows <= A_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(work_size);

    assert(success);

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(work_size);

        std::vector<LAPACK_int> pivots_local;

        const std::size_t sizeof_index_type = sizeof(index_type);
        const std::size_t sizeof_LAPACK_int = sizeof(LAPACK_int);

        LAPACK_int* pivots_ptr;

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4127 ) // conditional expression is constant
#endif
        if(sizeof_LAPACK_int != sizeof_index_type)
#ifdef _MSC_VER
#pragma warning( pop )
#endif
        {
            try
            {
                pivots_local.resize(num_cols);
            }
            catch(const std::exception& exc)
            {
                assert(false);

                internal_api_error_set_last(
                    (std::string("dense_matrix_QR_pivoted: Exception. ") + exc.what()));

                return false;
            }

            pivots_ptr = &pivots_local.front();
        }
        else
        {
            pivots_ptr = reinterpret_cast<LAPACK_int*>(pivots);
        }

        int info = 0;
        LAPACK_geqp3(
            &LAPACK_m, &LAPACK_n, A_col_values, &LAPACK_lda,
            pivots_ptr, tau_reflectors,
            work, &LAPACK_wsz, rwork,
            &info);

        success = (info == 0);

        assert(success);

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4127 ) // conditional expression is constant
#endif
        if(success && sizeof_LAPACK_int != sizeof_index_type)
#ifdef _MSC_VER
#pragma warning( pop )
#endif
        {
            for(index_type j = 0; (j < num_cols) && success; ++j)
            {
                success = integral_type_range_check_val<index_type>::in_non_negative_range(pivots_ptr[j]);
                assert(success);
                pivots[j] = index_type(pivots_ptr[j]);
            }
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_QR_pivoted: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
index_type dense_matrix_QR_pivoted_small_count(
    index_type  num_rows,
    index_type  num_cols,
    value_type* pivoted_QR_col_values,  // pivoted QR output, not original "A"
    index_type  pivoted_QR_col_leading_dim,
    typename precision_traits<value_type>::scalar fuzz)
{
    bool success =
        pivoted_QR_col_values &&
        num_rows <= pivoted_QR_col_leading_dim &&
        0 < fuzz;

    index_type small_count = std::numeric_limits<index_type>::max();

    if(success)
    {
        typedef typename precision_traits<value_type>::scalar precision_scalar;

        const index_type min_row_col = std::min(num_rows, num_cols);

        const precision_scalar max_abs_diag = std::abs(*pivoted_QR_col_values);

        const precision_scalar threshold = fuzz * precision_scalar(min_row_col) *
            std::numeric_limits<precision_scalar>::epsilon() * max_abs_diag;

        small_count = std_utils_count_less_equal_abs_reverse_inc(
            pivoted_QR_col_values,
            min_row_col,
            std::size_t(pivoted_QR_col_leading_dim) + 1,
            threshold);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_QR_pivoted_small_count: Error.");

    return small_count;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
index_type dense_matrix_QR_pivoted_right_null_size(
    index_type  num_rows,
    index_type  num_cols,
    value_type* pivoted_QR_col_values,  // pivoted QR output, not original "A"
    index_type  pivoted_QR_col_leading_dim,
    typename precision_traits<value_type>::scalar fuzz)
{
    const index_type small_count = dense_matrix_QR_pivoted_small_count(
        num_rows, num_cols,
        pivoted_QR_col_values, pivoted_QR_col_leading_dim,
        fuzz);

    const index_type right_null_size =
        num_rows > num_cols ? small_count : index_type(small_count + (num_cols - num_rows));

    return right_null_size;
}

// -----------------------------------------------------------------------------

// Call this function to get lwork for GEQRF.

template<typename index_type, typename value_type>
std::size_t dense_matrix_QR_lwork(
    index_type num_rows,
    index_type num_cols,
    index_type A_col_leading_dim,
    value_type /**/)
{
    bool success =
        num_rows <= A_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim);

    assert(success);

    value_type work = value_type();

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(-1);

        int info = 0;
        LAPACK_geqrf(
            &LAPACK_m, &LAPACK_n,
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
            "dense_matrix_QR_lwork: Error.");

    return success ? std::size_t(std::real(work)) : std::numeric_limits<std::size_t>::max();
}

// -----------------------------------------------------------------------------

// This function exists solely to wrap LAPACK GEQRF and do some common checks.

template<typename index_type, typename value_type>
bool dense_matrix_QR(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type* tau_reflectors,
    value_type* work,
    std::size_t work_size)
{
    bool success =
        A_col_values &&
        (std::min(num_rows, num_cols) == 0 || tau_reflectors) &&
        work &&
        num_rows <= A_col_leading_dim &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_rows) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(num_cols) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(A_col_leading_dim) &&
        integral_type_range_check_val<LAPACK_int>::in_non_negative_range(work_size);

    assert(success);

    if(success)
    {
        LAPACK_int LAPACK_m   = LAPACK_int(num_rows);
        LAPACK_int LAPACK_n   = LAPACK_int(num_cols);
        LAPACK_int LAPACK_lda = LAPACK_int(A_col_leading_dim);
        LAPACK_int LAPACK_wsz = LAPACK_int(work_size);

        int info = 0;
        LAPACK_geqrf(
            &LAPACK_m, &LAPACK_n, A_col_values, &LAPACK_lda,
            tau_reflectors,
            work, &LAPACK_wsz,
            &info);

        success = (info == 0);
        assert(success);
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_QR: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Use this function only if it is known that the matrix has full column-rank.
// Returns the column-space using unpivoted QR.

template<typename index_type, typename value_type>
bool dense_matrix_QR_orth_col_space_for_full_rank(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim)
{
    const index_type num_reflectors = std::min(num_rows, num_cols);

    const std::size_t work_size1 = dense_matrix_QR_lwork<index_type, value_type>(
            num_rows,
            num_cols,
            A_col_leading_dim,
            value_type());

    const std::size_t work_size2 = dense_matrix_reflectors_to_orth_col_lwork<index_type, value_type>(
            num_rows,
            num_cols,
            num_reflectors,
            A_col_leading_dim,
            value_type());

    const std::size_t work_size = std::max(work_size1, work_size2);

    std::vector<value_type> tau_reflectors;
    std::vector<value_type> work;

    try
    {
        tau_reflectors.resize(num_reflectors);
        work.resize(work_size);
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("dense_matrix_QR_orth_col_space_for_full_rank: Exception. ") + exc.what()));

        return false;
    }

    bool success =
        dense_matrix_QR(
            num_rows,
            num_cols,
            A_col_values,
            A_col_leading_dim,
            tau_reflectors.size() ? &tau_reflectors.front() : reinterpret_cast<value_type*>(0),
            work.size() ? &work.front() : reinterpret_cast<value_type*>(0),
            work_size)
        &&
        dense_matrix_reflectors_to_orth_col(
            num_rows,
            num_cols,
            num_reflectors,
            A_col_values,
            A_col_leading_dim,
            tau_reflectors.size() ? &tau_reflectors.front() : reinterpret_cast<value_type*>(0),
            work.size() ? &work.front() : reinterpret_cast<value_type*>(0),
            work_size);

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_QR_orth_col_space_for_full_rank: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_QR_H
