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


#ifndef SPARSE_VECTORS_UNION_W_TRANS_H
#define SPARSE_VECTORS_UNION_W_TRANS_H

// -----------------------------------------------------------------------------

#include "sparsity_union/sparse_vectors_union.h"
#include "sparse_vectors/sparse_vectors_transpose.h"
#include "sparse_vectors/sparse_vectors.h"
#include "internal_api_error/internal_api_error.h"
#include <vector>
#include <cstddef>

// -----------------------------------------------------------------------------
// Objective: To compute union of a sparse id collection of vectors with the
// transpose of another given sparse id collection of vectors (or self).
// -----------------------------------------------------------------------------

// ids in each vec of ids_vecs must be sorted.

// Search for Template-concept in the source to see requirements of
// ids_collection_type and inv_id_func_collection_type.

template
<
    typename index_type,
    typename ids_collection_type,
    typename inv_id_func_collection_type
>
bool sparse_vectors_union_w_trans(
    const ids_collection_type& vecs,
    const inv_id_func_collection_type& ids_vecs,
    const ids_collection_type& vecs_for_trans,
    std::vector< std::vector<index_type> >& trans_vecs_union)
{
    // Choose largest possible offset_size (std::size_t).
    sparse_vectors_ids<index_type, std::size_t> trans_of_vecs_for_trans;

    bool success =
        sparse_vectors_transpose_ids(
            vecs_for_trans.num_vecs(),
            vecs_for_trans,
            ids_vecs,   // trans of trans is identity
            trans_of_vecs_for_trans)
            &&
        sparse_vectors_union(
            vecs_for_trans.num_vecs(),
            vecs,
            trans_of_vecs_for_trans,
            trans_vecs_union);

   assert(success);

    if(!success)
        internal_api_error_set_last("sparse_vectors_union_w_trans: Error.");

    return success;
}

// -----------------------------------------------------------------------------

// ids in each vec of ids_vecs must be sorted.

template
<
    typename index_type,
    typename ids_collection_type,
    typename inv_id_func_collection_type
>
bool sparse_vectors_union_w_self_trans(
    const ids_collection_type& vecs,
    const inv_id_func_collection_type& ids_vecs,
    std::vector< std::vector<index_type> >& trans_vecs_union)
{
    bool success = sparse_vectors_union_w_trans(
        vecs,
        ids_vecs,
        vecs,
        trans_vecs_union);

   assert(success);

    if(!success)
        internal_api_error_set_last("sparse_vectors_union_w_self_trans: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_VECTORS_UNION_W_TRANS_H
