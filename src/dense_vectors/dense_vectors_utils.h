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


#ifndef DENSE_VECTORS_UTILS_H
#define DENSE_VECTORS_UTILS_H

// -----------------------------------------------------------------------------

#include "math/vector_utils.h"
#include "internal_api_error/internal_api_error.h"
#include <algorithm>   // std::{fill, copy}
#include <cstddef>
#include <cassert>

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_vectors_utils_axpby(
    index_type  num_vecs,
    index_type  max_size,
    const value_type* X_values,
    index_type  X_leading_dim,
    value_type* Y_values,
    index_type  Y_leading_dim,
    value_type  a,
    value_type  b)
{
    // Y <- a * X + b * Y

    bool success =
        X_values &&
        Y_values &&
        max_size <= X_leading_dim &&
        max_size <= Y_leading_dim;

   assert(success);

    if(success)
    {
        const value_type* X_j = X_values;
        value_type* Y_j = Y_values;

        for(index_type j = 0; j < num_vecs;
            ++j,
            X_j += X_leading_dim,
            Y_j += Y_leading_dim)
        {
            vector_utils_axpby(max_size, X_j, Y_j, a, b,
                index_type(1), index_type(1));
        }
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("dense_vectors_utils_axpby: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_vectors_utils_fill(
    index_type  num_vecs,
    index_type  max_size,
    value_type* A_values,
    index_type  A_leading_dim,
    value_type  fill_with)
{
    bool success =
        A_values &&
        max_size <= A_leading_dim;

   assert(success);

    if(success)
    {
        const value_type* A_j_end = A_values +
            std::size_t(A_leading_dim) * std::size_t(num_vecs);

        for(value_type* A_j = A_values;
            A_j < A_j_end;
            A_j += A_leading_dim)
        {
            std::fill(A_j, A_j + max_size, fill_with);
        }
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("dense_vectors_utils_fill: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool dense_vectors_utils_copy(
    index_type num_vecs,
    index_type max_size,
    const value_type* A_values,
    index_type A_leading_dim,
    value_type* B_values,
    index_type B_leading_dim)
{
    bool success =
        max_size <= A_leading_dim &&
        max_size <= B_leading_dim &&
        A_values &&
        B_values;

   assert(success);

    if(success)
    {
        const value_type* values_it = A_values;
        value_type* B_it = B_values;

        for(index_type i = 0; i < num_vecs; ++i,
            values_it += A_leading_dim,
            B_it += B_leading_dim)
        {
            std::copy(values_it, values_it + max_size, B_it);
        }
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("dense_vectors_utils_copy: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // DENSE_VECTORS_UTILS_H
