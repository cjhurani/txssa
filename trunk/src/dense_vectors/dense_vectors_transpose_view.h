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


#ifndef DENSE_VECTORS_TRANSPOSE_VIEW_H
#define DENSE_VECTORS_TRANSPOSE_VIEW_H

// -----------------------------------------------------------------------------

#include "dense_vectors/dense_vectors.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Define simple C++ adaptor class to wrap transpose of dense vectors
// without any allocation.
// -----------------------------------------------------------------------------

// Transposed view of dense vectors with uniform offsets (between their
// "begin").  No data allocation happens for the transpose.

template<typename index_t, typename value_t>
class dense_vectors_transpose_view
{
public:

    typedef index_t index_type;
    typedef value_t value_type;

    // id_func_type and inv_id_func_type are the same.
    typedef dense_vectors_id_func<index_type> id_func_type;
    typedef dense_vectors_id_func<index_type> inv_id_func_type;

    dense_vectors_transpose_view(
        const dense_vectors<index_type, value_type>& in_dv)
        :
        n_vecs(in_dv.num_vecs()),
        each_vec_size(in_dv.vec_size()),
        leading_dim(in_dv.leading_dimension()),
        vec_values(in_dv.vec_values())
    {
        assert(vec_values);
        assert(0 < leading_dim);
        assert(each_vec_size <= leading_dim);
    }

    index_type num_vecs() const
    {
        return each_vec_size; // flip
    }

    index_type vec_size() const
    {
        return n_vecs; // flip
    }

    // Compatibility with sparse vectors
    index_type max_size() const
    {
        return n_vecs; // flip
    }

    // Compatibility with sparse vectors
    index_type num_vec_entries(index_type i) const
    {
        (void) i;
        assert(i < each_vec_size); // flip

        return n_vecs; // flip
    }

    value_type* vec_values_begin(index_type i)
    {
        assert(i < each_vec_size); // flip

        return vec_values + i;
    }

    const value_type* vec_values_begin(index_type i) const
    {
        assert(i < each_vec_size); // flip

        return vec_values + i;
    }

    id_func_type id_func(index_type i) const
    {
        (void) i;
        assert(i < each_vec_size); // flip

        return id_func_type(n_vecs); // flip
    }

    inv_id_func_type inv_id_func(index_type i) const
    {
        (void) i;
        assert(i < each_vec_size); // flip

        return inv_id_func_type(n_vecs); // flip
    }

    index_type inc(index_type i) const
    {
        (void) i;
        // Note: this is inc in memory, not in mathematical component id.

        assert(i < each_vec_size); // flip

        return leading_dim;
    }

private:

    index_type  n_vecs;
    index_type  each_vec_size;
    index_type  leading_dim;
    value_type* vec_values;
};

// -----------------------------------------------------------------------------

#endif // DENSE_VECTORS_TRANSPOSE_VIEW_H
