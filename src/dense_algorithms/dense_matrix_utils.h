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


#ifndef DENSE_MATRIX_UTILS_H
#define DENSE_MATRIX_UTILS_H

// -----------------------------------------------------------------------------

#include "dense_vectors/dense_vectors_utils.h"
#include "math/complex_types.h"
#include "cpp/std_extensions.h"
#include "internal_api_error/internal_api_error.h"
#include <algorithm>   // std::{fill, min, copy}
#include <complex>
#include <cassert>

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_diagonal_add(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type  to_add)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < min_row_col; ++j)
        {
            A_j[j] += to_add;
            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_diagonal_add: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_fill_strict_lower(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type  fill_with)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < min_row_col; ++j)
        {
            std::fill(A_j + j + 1, A_j + num_rows, fill_with);
            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_fill_strict_lower: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_fill_strict_upper(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type  fill_with)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 1; j < min_row_col; ++j)
        {
            std::fill(A_j, A_j + j, fill_with);
            A_j += A_col_leading_dim;
        }

        // Fill everything if more columns left.
        if(num_rows < num_cols)
        {
            success = dense_vectors_utils_fill(
                index_type(num_cols - num_rows),
                num_rows,
                A_j,
                A_col_leading_dim,
                fill_with);

            assert(success);
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_fill_strict_upper: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_fill_lower(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type  fill_with)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < min_row_col; ++j)
        {
            std::fill(A_j + j, A_j + num_rows, fill_with);
            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_fill_lower: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_fill_upper(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type  fill_with)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim;

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < min_row_col; ++j)
        {
            std::fill(A_j, A_j + j + 1, fill_with);
            A_j += A_col_leading_dim;
        }

        // Fill everything if more columns left.
        if(num_rows < num_cols)
        {
            success = dense_vectors_utils_fill(
                index_type(num_cols - num_rows),
                num_rows,
                A_j,
                A_col_leading_dim,
                fill_with);

            assert(success);
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_fill_upper: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<    typename index_type,
    typename value_type,
    typename unary_func
>
bool dense_matrix_utils_copy_lower_to_upper_func(
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    const unary_func& func)
{
    bool success =
        A_col_values &&
        matrix_size <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < matrix_size; ++j)
        {
            for(index_type i = j + 1; i < matrix_size; ++i)
            {
                // (i,j) is in strict lower half.
                // A(j,i) = func(A(i,j)) for j < i.
                A_col_values[std::size_t(A_col_leading_dim) * std::size_t(i) + std::size_t(j)] = func(A_j[i]);
            }

            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_copy_lower_to_upper_func: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<    typename index_type,
    typename value_type,
    typename unary_func
>
bool dense_matrix_utils_copy_upper_to_lower_func(
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    const unary_func& func)
{
    bool success =
        A_col_values &&
        matrix_size <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        const value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < matrix_size; ++j)
        {
            for(index_type i = 0; i < j; ++i)
            {
                // (i,j) is in strict upper half.
                // A(j,i) = func(A(i,j)) for i < j.
                A_col_values[std::size_t(A_col_leading_dim) * std::size_t(i) + std::size_t(j)] = func(A_j[i]);
            }

            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_copy_upper_to_lower_func: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_transpose_in_place(
    index_type  matrix_size,
    value_type* A_col_values,
    index_type  A_col_leading_dim)
{
    bool success =
        A_col_values &&
        matrix_size <= A_col_leading_dim;

    assert(success);

    if(success)
    {
        value_type* A_j = A_col_values; // beginning of jth column.

        for(index_type j = 0; j < matrix_size; ++j)
        {
            for(index_type i = 0; i < j; ++i)
            {
                // (i,j) is in strict upper half.
                const std::size_t lower_id = std::size_t(A_col_leading_dim) * std::size_t(i) + std::size_t(j);
                const value_type upper_conj = std::conj(A_j[i]);
                A_j[i] = std::conj(A_col_values[lower_id]);
                A_col_values[lower_id] = upper_conj;
            }

            // For complex, also need to take conj of diagonal value.

            A_j[j] = std::conj(A_j[j]);

            A_j += A_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_transpose_in_place: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_matrix_utils_copy_strict_lower(
    index_type  num_rows,
    index_type  num_cols,
    const value_type* A_col_values,
    index_type  A_col_leading_dim,
    value_type* B_col_values,
    index_type  B_col_leading_dim)
{
    bool success =
        A_col_values &&
        B_col_values &&
        num_rows <= A_col_leading_dim &&
        num_rows <= B_col_leading_dim;

    assert(success);

    if(success)
    {
        const index_type min_row_col = std::min(num_rows, num_cols);
        const value_type* A_j = A_col_values; // beginning of jth column.
        value_type* B_j = B_col_values;       // beginning of jth column.

        for(index_type j = 0; j < min_row_col; ++j)
        {
            std::copy(A_j + j + 1, A_j + num_rows, B_j + j + 1);

            A_j += A_col_leading_dim;
            B_j += B_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_copy_strict_lower: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// Assume that ATA is generated from A that is complex-symmetric
template<typename index_type, typename scalar_type>
bool dense_matrix_utils_complex_sym_compute_AAT_from_ATA(
    index_type size,
    const std::complex<scalar_type>* ATA_col_values,
    index_type ATA_col_leading_dim,
    std::complex<scalar_type>* AAT_col_values,
    index_type AAT_col_leading_dim)
{
    bool success =
        ATA_col_values &&
        AAT_col_values &&
        size <= ATA_col_leading_dim &&
        size <= AAT_col_leading_dim;

    // ATA == conj(AAT)

    assert(success);

    if(success)
    {
        for(index_type j = 0; j < size; ++j)
        {
            for(index_type i = 0; i < size; ++i)
            {
                AAT_col_values[i] = std::conj(ATA_col_values[i]);
            }

            ATA_col_values += ATA_col_leading_dim;
            AAT_col_values += AAT_col_leading_dim;
        }
    }

    if(!success)
        internal_api_error_set_last(
            "dense_matrix_utils_complex_sym_compute_AAT_from_ATA: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_MATRIX_UTILS_H
