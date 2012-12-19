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


#ifndef P_NORM_SPARSITY_VECTORS_H
#define P_NORM_SPARSITY_VECTORS_H

// -----------------------------------------------------------------------------

#include "p_norm_sparsity_vectors/p_norm_sparsity_vector.h"
#include "internal_api_error/internal_api_error.h"
#include <string>
#include <vector> 
#include <algorithm>   // std::copy
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: The functions below compute the sparsity patterns for a collection
// of vectors.
//
// The patterns for a matrix require combining row and column patterns after
// finding patterns for both.  Those operations are not done here.  Hence, the
// concept of working with collection of vectors is used here.  By choosing
// appropriate vals_inc_collection_type, one can make this function work for
// almost any matrix - row/col-oriented dense or CSR, or CSC.
//
// For sparse vectors, the offset ids are stored in the pattern and not absolute
// position of a non-zero in the vector.
//
// The reason we return the sparsity pattern using a vector of vectors and not
// something like "CSR" format is because the length of each sparse vectors
// (and the global array) is not known a-priori, so dealing with vector of
// vectors is simpler.
// -----------------------------------------------------------------------------

// Template-concept: vals_inc_collection_type
// Following member functions should be present.
//
// index_type          num_vecs        ()           const
// [const] value_type* vec_values_begin(index_type) const
// index_type          num_vec_entries (index_type) const
// index_type          inc             (index_type) const

template
<
    typename index_type,
    typename scalar_type,
    typename vals_inc_collection_type
>
bool p_norm_sparsity_vectors(

// algorithmic options:
    scalar_type ratio,
    scalar_type p,
    index_type min_num_nnz,

// input data:
    index_type max_vec_size,
    const vals_inc_collection_type& vecs,

// output ids:
    std::vector< std::vector<index_type> >& vec_ids)
{
    typedef typename vals_inc_collection_type::value_type value_type;

    const index_type num_vecs = vecs.num_vecs();

    std::vector<index_type> out_ids;
    std::vector<scalar_type> work_val;
    std::vector< std::vector<index_type> > tmp_vec_ids;

    bool success = false;

    try
    {
        // Temporary allocation, could throw.
        out_ids.resize(max_vec_size);
        work_val.resize(max_vec_size);
        tmp_vec_ids.resize(num_vecs);

        success = true;

        for(index_type i = 0; i < num_vecs; ++i)
        {
            const index_type num_vec_entries = vecs.num_vec_entries(i);
            const value_type* values_begin = vecs.vec_values_begin(i);

            assert(num_vec_entries <= max_vec_size);
            assert(values_begin);

            index_type out_num_nnz;
            success = p_norm_sparsity_vector(
                ratio, p, min_num_nnz,
                num_vec_entries, values_begin, vecs.inc(i),
                &out_num_nnz, &out_ids.front(),
                &work_val.front());

            if(success)
            {
                tmp_vec_ids[i].resize(out_num_nnz); // Could throw
                std::copy(out_ids.begin(), out_ids.begin() + out_num_nnz,
                    tmp_vec_ids[i].begin());
            }
            else
            {
                break;
            }
        }

        if(success)
            vec_ids.swap(tmp_vec_ids);
    }
    catch(const std::exception& exc)
    {
        internal_api_error_set_last(
            (std::string("p_norm_sparsity_vectors: Exception. ") + exc.what()));

        assert(false);
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last("p_norm_sparsity_vectors: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_VECTORS_H
