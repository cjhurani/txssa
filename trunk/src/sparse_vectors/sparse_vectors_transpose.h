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


#ifndef SPARSE_VECTORS_TRANSPOSE_H
#define SPARSE_VECTORS_TRANSPOSE_H

// -----------------------------------------------------------------------------

#include "internal_api_error/internal_api_error.h"
#include <vector>
#include <string>
#include <algorithm>   // std::fill
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Transpose sparsity pattern of a collection of vectors.  The ids
// could be absolute ids or relative ids.  Hence one also requires a collection
// of "parent" ids to be available.
// -----------------------------------------------------------------------------

// Template-concept: ids_collection_type
// Following member functions should be present.
//
// index_type          num_vecs        ()           const
// [const] index_type* vec_ids_begin   (index_type) const
// index_type          num_vec_entries (index_type) const
//
// Template-concept: inv_id_func_collection_type
// Following member functions should be present.
//
// index_type num_vecs() const
// inv_id_func_collection_type::inv_id_func_type inv_id_func(index_type) const
//
// Template-concept: out_collection_type
// Following member functions should be present.
//
// Default constructor
// Destructor
// bool        allocate     (index_type, index_type, [const] index_type*)
// index_type* vec_ids_begin(index_type)
// [void]      swap(out_collection_type&)
//
// Template-concept: rand_it
// Following should be possible (or exist) for variables of type rand_it.
// Copy-able
// index_type& operator[](index_type)
//
// Template-concept: id_func inv_id_func should support
// index_type operator()(index_type) const
// either as a functor or global functions.

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename ids_collection_type,
    typename rand_it
>
void sparse_vectors_transpose_nnz_add(
    index_type num_vecs,
    const ids_collection_type& vecs,
    rand_it transpose_nnz)
{
    assert(num_vecs == vecs.num_vecs());

    for(index_type i = 0; i < num_vecs; ++i)
    {
        const index_type num_vec_entries = vecs.num_vec_entries(i);
        const index_type* vec            = vecs.vec_ids_begin(i);
        assert(num_vec_entries <= vecs.num_vec_entries(i));
        assert(vec);

        for(index_type j = 0; j < num_vec_entries; ++j)
        {
            ++transpose_nnz[vec[j]];
        }
    }
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename ids_collection_type,         // What we have to transpose
    typename inv_id_func_collection_type, // To convert from "actual" in trans
    typename out_collection_type
>
bool sparse_vectors_transpose_ids(
    index_type num_vecs,
    const ids_collection_type& vecs,
    const inv_id_func_collection_type& trans_ids_vecs,
    out_collection_type& trans_vecs)
{
    if(num_vecs != vecs.num_vecs())
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_vectors_transpose_ids: num_vecs != vecs.num_vecs().");

        return false;
    }

    const index_type max_vec_size = trans_ids_vecs.num_vecs();

    std::vector<index_type> transpose_nnz;

    out_collection_type tmp_trans_vecs;

    try
    {
        transpose_nnz.resize(max_vec_size);

        sparse_vectors_transpose_nnz_add(
            num_vecs, vecs, transpose_nnz.begin());

        bool success =
            tmp_trans_vecs.allocate(max_vec_size, num_vecs, &transpose_nnz.front());

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors_transpose_ids: Allocation error.");

            return false;
        }
    }
    catch(const std::exception& exc)
    {
        assert(false);

        internal_api_error_set_last(
            (std::string("sparse_vectors_transpose_ids: Exception. ") + exc.what()));

        return false;
    }

    // Reuse space.  "Rename" first.
    std::vector<index_type>& ids_done = transpose_nnz;
    std::fill(ids_done.begin(), ids_done.end(), index_type(0));

    for(index_type i = 0; i < num_vecs; ++i)
    {
        const index_type num_vec_entries = vecs.num_vec_entries(i);
        const index_type* ids = vecs.vec_ids_begin(i);
        assert(num_vec_entries <= max_vec_size);
        assert(ids);

        for(index_type j = 0; j < num_vec_entries; ++j)
        {
            const index_type id_mapped = ids[j];

            // Now we cannot put "i" in tmp_trans_vecs. This is because "i" is
            // component id and not necessarily derived from trans_ids_vecs.
            // Find i in id_mapped vector of trans_ids_vecs and use that.

            const index_type i_off = trans_ids_vecs.inv_id_func(id_mapped)(i);

            tmp_trans_vecs.vec_ids_begin(id_mapped)[ids_done[id_mapped]++] =
                i_off;
        }
    }

    trans_vecs.swap(tmp_trans_vecs);

    return true;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_VECTORS_TRANSPOSE_H
