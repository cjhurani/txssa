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


#ifndef SPARSE_VECTORS_UNION_H
#define SPARSE_VECTORS_UNION_H

// -----------------------------------------------------------------------------

#include "cpp/std_new_features.h"     // std_new_features_is_sorted
#include "internal_api_error/internal_api_error.h"
#include <string>
#include <vector>
#include <algorithm>                  // std::{copy, sort, set_union}
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: To compute union of a sparse id collection of vectors with another
// given sparse id collection of vectors.
// -----------------------------------------------------------------------------

// Template-concept: ids_collection_type_1 and _2 model ids_collection_type.

template
<
    typename index_type,
    typename ids_collection_type_1,
    typename ids_collection_type_2
>
bool sparse_vectors_union(
    index_type max_vec_size,
    const ids_collection_type_1& vecs_1,
    const ids_collection_type_2& vecs_2,
    std::vector< std::vector<index_type> >& vecs_union)
{
    const index_type num_vecs = vecs_1.num_vecs();

    bool success = (vecs_1.num_vecs() == vecs_2.num_vecs());

    if(!success)
    {
        assert(false);

        internal_api_error_set_last("sparse_vectors_union: vecs_1.num_vecs() != vecs_2.num_vecs().");

        return success;
    }

    success = false;

    typedef typename std::vector<index_type>::size_type size_type;

    try
    {
        std::vector< std::vector<index_type> > tmp_vecs_union(num_vecs);
    
        std::vector<index_type>
            tmp_1(max_vec_size),
            tmp_2(max_vec_size),
            tmp_12(2 * size_type(max_vec_size));

        index_type* vec_1_beg = &tmp_1.front();
        index_type* vec_2_beg = &tmp_2.front();

        for(index_type i = 0; i < num_vecs; ++i)
        {
            const index_type num_vec_entries_1 = vecs_1.num_vec_entries(i);
            const index_type num_vec_entries_2 = vecs_2.num_vec_entries(i);

            const index_type* vec_1 = vecs_1.vec_ids_begin(i);
            const index_type* vec_2 = vecs_2.vec_ids_begin(i);

            assert(num_vec_entries_1 <= max_vec_size);
            assert(num_vec_entries_2 <= max_vec_size);
            assert(vec_1);
            assert(vec_2);

            const bool sorted_1 = std_new_features_is_sorted(
                vec_1, vec_1 + num_vec_entries_1);

            if(!sorted_1)
            {
                std::copy(vec_1, vec_1 + num_vec_entries_1, vec_1_beg);
                std::sort(vec_1_beg, vec_1_beg + num_vec_entries_1);
            }

            const index_type* vec_1_beg_const = sorted_1 ? vec_1 : vec_1_beg;

            const bool sorted_2 = std_new_features_is_sorted(
                vec_2, vec_2 + num_vec_entries_2);

            if(!sorted_2)
            {
                std::copy(vec_2, vec_2 + num_vec_entries_2, vec_2_beg);
                std::sort(vec_2_beg, vec_2_beg + num_vec_entries_2);
            }

            const index_type* vec_2_beg_const = sorted_2 ? vec_2 : vec_2_beg;

            const size_type union_size =
                std::set_union(
                    vec_1_beg_const, vec_1_beg_const + num_vec_entries_1,
                    vec_2_beg_const, vec_2_beg_const + num_vec_entries_2,
                    tmp_12.begin()) -
                    tmp_12.begin();

            tmp_vecs_union[i].resize(union_size); // Could throw

            std::copy(tmp_12.begin(), tmp_12.begin() + union_size,
                tmp_vecs_union[i].begin());
        }

        vecs_union.swap(tmp_vecs_union);

        success = true;
    }
    catch(const std::exception& exc)
    {
        internal_api_error_set_last(
            (std::string("sparse_vectors_union: Exception. ") + exc.what()));

        assert(false);
    }

   assert(success);

    if(!success)
        internal_api_error_set_last("sparse_vectors_union: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_VECTORS_UNION_H
