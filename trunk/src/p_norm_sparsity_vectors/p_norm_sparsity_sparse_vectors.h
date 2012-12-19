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


#ifndef P_NORM_SPARSITY_SPARSE_VECTORS_H
#define P_NORM_SPARSITY_SPARSE_VECTORS_H

// -----------------------------------------------------------------------------

#include "p_norm_sparsity_vectors/p_norm_sparsity_vectors.h"
#include "sparse_vectors/sparse_vectors.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <sstream>
#include <vector>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: The functions below compute sparsity pattern (as a vector of
// patterns) for collections of sparse vectors.
// -----------------------------------------------------------------------------

// The ids in output vectors are not necessarily sorted.  Besides, they are
// relative ids.  To use them to access vector value_type, you need the
// underlying offset ids.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool p_norm_sparsity_sparse_vectors(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio,
    typename precision_traits<value_type>::scalar p,
    index_type min_num_nnz,

// input data:
    index_type num_vecs,
    index_type max_vec_size,
    const offset_type* vec_offsets,
    const value_type* vec_values,

// output ids:
    std::vector< std::vector<index_type> >& vec_ids)
{
    if(!vec_offsets || !vec_values)
    {
        assert(false);

        std::ostringstream oss;
        oss << "p_norm_sparsity_sparse_vectors: Unacceptable input argument(s).";
        if(!vec_offsets) oss << " !vec_offsets.";
        if(!vec_values)  oss << " !vec_values.";

        internal_api_error_set_last(oss.str());

        return false;
    }

    sparse_vectors<const index_type, const offset_type, const value_type> vecs(
        num_vecs,
        max_vec_size,
        vec_offsets,
        (const index_type*) 0, // Note: no need for ids
        vec_values);

    bool success = p_norm_sparsity_vectors(
        ratio,
        p,
        min_num_nnz,
        max_vec_size,
        vecs,
        vec_ids);

   assert(success);

    if(!success)
        internal_api_error_set_last("p_norm_sparsity_sparse_vectors: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_SPARSE_VECTORS_H
