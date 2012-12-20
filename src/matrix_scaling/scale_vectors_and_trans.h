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


#ifndef SCALE_VECTORS_AND_TRANS_H
#define SCALE_VECTORS_AND_TRANS_H

// -----------------------------------------------------------------------------

#include <cassert>     // assert

// -----------------------------------------------------------------------------
// Objective: Multiply entries in a collection of vectors by scales depending
// on vector id and entry's position id.  The name "scale_vectors_and_trans"
// means one is scaling each vector by a single number as well as scaling the
// transposes by a single number.
// -----------------------------------------------------------------------------

// Template-concept: id_func_collection_type
// Following member functions should be present.
//
// index_type                        num_vecs       ()           const
// index_type                        num_vec_entries(index_type) const
// id_func_collection_type::id_func_type id_func    (index_type) const

// Template-concept: vals_id_func_collection_type
// Following member functions should be present.
//
// index_type                            num_vecs        ()           const
// index_type                            num_vec_entries (index_type) const
// id_func_collection_type::id_func_type id_func         (index_type) const
// value_type*                           vec_values_begin(index_type) const

template
<
    typename index_type,
    typename scale_type,
    typename vals_id_func_collection_type
>
void scale_vectors_and_trans(
    const scale_type* scales,
    const scale_type* trans_scales,
    vals_id_func_collection_type& vecs)
{
    assert(scales);
    assert(trans_scales);

    const index_type num_vecs = vecs.num_vecs();

    typedef typename vals_id_func_collection_type::id_func_type id_func_type;
    typedef typename vals_id_func_collection_type::value_type   value_type;

    for(index_type i = 0; i < num_vecs; ++i)
    {
        const index_type num_vec_entries = vecs.num_vec_entries(i);
        value_type* vec_vals             = vecs.vec_values_begin(i);

        const id_func_type& id_func = vecs.id_func(i);

        const scale_type i_scale = scales[i];

        for(index_type j = 0; j < num_vec_entries; ++j)
        {
            vec_vals[j] *= i_scale * trans_scales[id_func(j)];
        }
    }
}

// -----------------------------------------------------------------------------

// Scale a collection with the same scale and trans_scale.  The collection must
// be square.

template
<
    typename index_type,
    typename scale_type,
    typename vals_id_func_collection_type
>
void scale_vectors_and_trans(
    const scale_type* scales,
    vals_id_func_collection_type& vecs)
{
    scale_vectors_and_trans
        <index_type, scale_type, vals_id_func_collection_type>
            (scales, scales, vecs);
}

// -----------------------------------------------------------------------------

#endif // SCALE_VECTORS_AND_TRANS_H
