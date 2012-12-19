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


#ifndef P_NORM_SPARSITY_DENSE_VECTORS_H
#define P_NORM_SPARSITY_DENSE_VECTORS_H

// -----------------------------------------------------------------------------

#include "p_norm_sparsity_vectors/p_norm_sparsity_vectors.h"
#include "dense_vectors/dense_vectors_transpose_view.h"
#include "dense_vectors/dense_vectors.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <sstream>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: The functions below compute sparsity pattern (as a vector of
// patterns) for collections of dense vectors.
// -----------------------------------------------------------------------------

// The ids in output vectors are not necessarily sorted.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_vectors(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz,

// input data:
    index_type num_vecs,
    index_type vec_size,
    index_type leading_dim,
    const value_type* vec_values,

// output ids:
    std::vector< std::vector<index_type> >& vec_ids)
{
    if(!vec_values || leading_dim < vec_size)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_vectors: Unacceptable input argument(s).";
        if(!vec_values)            oss << " !vec_values.";
        if(leading_dim < vec_size) oss << " leading_dim < vec_size.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    const dense_vectors<index_type, const value_type> dense_vecs(
        num_vecs,
        vec_size,
        leading_dim,
        vec_values);

    bool success = p_norm_sparsity_vectors(
        ratio,
        p,
        min_num_nnz,
        vec_size,
        dense_vecs,
        vec_ids);

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_vectors: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// The data is of original collection of dense vectors, but the pattern is for
// transposed set of vectors.

template<typename index_type, typename value_type>
bool p_norm_sparsity_dense_vectors_transpose_view(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz,

// input data, not transposed:
    index_type num_vecs,
    index_type vec_size,
    index_type leading_dim,
    const value_type* vec_values,

// output ids:
    std::vector< std::vector<index_type> >& vec_ids)
{
    if(!vec_values || leading_dim < vec_size)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_dense_vectors_transpose_view: Unacceptable input argument(s).";
        if(!vec_values)            oss << " !vec_values.";
        if(leading_dim < vec_size) oss << " leading_dim < vec_size.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    const dense_vectors<index_type, const value_type> dv(
        num_vecs,
        vec_size,
        leading_dim,
        vec_values);

    const dense_vectors_transpose_view<index_type, const value_type>
        transpose_view(dv);

    bool success = p_norm_sparsity_vectors(
        ratio,
        p,
        min_num_nnz,
        num_vecs,  // flip
        transpose_view,
        vec_ids);

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_dense_vectors_transpose_view: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_DENSE_VECTORS_H
