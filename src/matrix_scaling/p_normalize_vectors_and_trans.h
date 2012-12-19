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


#ifndef P_NORMALIZE_VECTORS_AND_TRANS_H
#define P_NORMALIZE_VECTORS_AND_TRANS_H

// -----------------------------------------------------------------------------

#include "scale_vectors_and_trans.h"
#include "p_norm_of_vectors/p_norm_of_vectors_and_trans.h"
#include "p_norm_of_vectors/p_norm_of_vectors.h"
#include "math/precision_traits.h"
#include "math/vector_utils.h"

#include <vector>
#include <cstddef>
#include <limits>
#include <algorithm>   // std::{transform, fill}
#include <functional>  // std::{ptr_fun, multiplies}
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Iteratively normalize collection of vectors so that each vec
// and each vec in transpose has p-norm equal to 1.
// -----------------------------------------------------------------------------

// Template-concept: id_func_collection_type

// Template-concept: vals_id_func_max_collection_type
// Following member functions should be present.
//
// index_type                            num_vecs        ()           const
// index_type                            max_size        ()           const
// index_type                            num_vec_entries (index_type) const
// id_func_collection_type::id_func_type id_func         (index_type) const
// value_type*                           vec_values_begin(index_type) const

template
<
    typename index_type,
    typename value_type,
    typename vals_id_func_max_collection_type
>
bool p_normalize_vectors_and_trans(

// algorithmic options:
    typename precision_traits<value_type>::scalar p,
    typename precision_traits<value_type>::scalar tolerance,
    std::size_t max_iters,

// input/output data:
    vals_id_func_max_collection_type& vecs,

// output data:
    typename precision_traits<value_type>::scalar* vec_scale,
    typename precision_traits<value_type>::scalar* trans_vec_scale,
    std::size_t* iters)
{
    typedef typename precision_traits<value_type>::scalar scalar_type;

    if(p < 1)
    {
        assert(false);
        return false;
    }

    const index_type num_vecs = vecs.num_vecs();
    const index_type max_size = vecs.max_size();

    std::size_t iters_taken = 0;
        
    if(0 < num_vecs && 0 < max_size)
    {
        index_type num_vecs_with_zero_norm = 0;
        index_type num_trans_vecs_with_zero_norm = 0;

        std::fill(      vec_scale,       vec_scale + num_vecs, scalar_type(1));
        std::fill(trans_vec_scale, trans_vec_scale + max_size, scalar_type(1));

        const scalar_type inf = std::numeric_limits<scalar_type>::infinity();
        const scalar_type zero = 0;

        scalar_type vec_limit, trans_vec_limit;

        if(p == inf)
        {
            vec_limit = 1;
        }
        else
        {
            vec_limit = std::pow(scalar_type(max_size)/num_vecs, scalar_type(0.25)/p);
        }

        trans_vec_limit = 1/vec_limit;

        std::vector<scalar_type> vec_norms, t_vec_norms;

        try
        {
            vec_norms.resize(num_vecs);
            t_vec_norms.resize(max_size);
        }
        catch(const std::exception& /*exc*/)
        {
            assert(false);
            return false;
        }

        scalar_type* vec_norms_ptr   = &vec_norms.front();
        scalar_type* t_vec_norms_ptr = &t_vec_norms.front();

        while(iters_taken < max_iters)
        {
            // First compute the norms and trans-norms.

            p_norm_of_vectors_and_trans
                <scalar_type, vals_id_func_max_collection_type>
                    (p, vecs, vec_norms_ptr, t_vec_norms_ptr);

            // We don't expect there will be any (or many) vectors or
            // trans-vectors whose norms will be zeros.  But if there are, we
            // want the algorithm to work without trouble.  Of course, the
            // counts and locations of such vectors will not change as
            // iterations proceed, so we do it once only.

            if(iters_taken == 0)
            {
                num_vecs_with_zero_norm = static_cast<index_type>(std::count(
                    vec_norms.begin(), vec_norms.end(), zero));

                num_trans_vecs_with_zero_norm = static_cast<index_type>(std::count(
                    t_vec_norms.begin(), t_vec_norms.end(), zero));

                if( num_vecs_with_zero_norm == num_vecs ||
                    num_trans_vecs_with_zero_norm == max_size)
                {
                    // If any one is true, the other must be true also.
                    assert(num_trans_vecs_with_zero_norm == max_size);
                    assert(num_vecs_with_zero_norm == num_vecs);

                    // all zero matrix, so break.
                    break;
                }

                // If we encountered zero norms, we have to change the limiting
                // values (if they depend on matrix size, which is when
                // p < inf).

                if( p != inf &&
                    (num_vecs_with_zero_norm > 0 ||
                     num_trans_vecs_with_zero_norm > 0))
                {
                    vec_limit = std::pow(
                        scalar_type(max_size - num_trans_vecs_with_zero_norm) /
                                  (num_vecs - num_vecs_with_zero_norm), scalar_type(0.25)/p);

                    trans_vec_limit = 1/vec_limit;
                }
            }

            std::transform(
                vec_norms.begin(), vec_norms.end(), vec_norms.begin(),
                std::ptr_fun<scalar_type, scalar_type>(std::sqrt));

            std::transform(
                t_vec_norms.begin(), t_vec_norms.end(), t_vec_norms.begin(),
                std::ptr_fun<scalar_type, scalar_type>(std::sqrt));

            if(tolerance >= 0 &&
                vector_utils_max_abs_diff_with_ignore(num_vecs, vec_norms_ptr,
                    vec_limit, zero) <= tolerance &&
                vector_utils_max_abs_diff_with_ignore(max_size, t_vec_norms_ptr,
                    trans_vec_limit, zero) <= tolerance)
            {
                break;
            }

            vector_utils_invert_non_zero(num_vecs,   vec_norms_ptr);
            vector_utils_invert_non_zero(max_size, t_vec_norms_ptr);

            scale_vectors_and_trans
                <index_type, scalar_type, vals_id_func_max_collection_type>
                    (vec_norms_ptr, t_vec_norms_ptr, vecs);

            std::transform(vec_norms.begin(), vec_norms.end(), vec_scale, vec_scale,
                std::multiplies<scalar_type>());

            std::transform(t_vec_norms.begin(), t_vec_norms.end(), trans_vec_scale, trans_vec_scale,
                std::multiplies<scalar_type>());

            ++iters_taken;
        }
    }

    if(iters)
        *iters = iters_taken;

    return true;
}

// -----------------------------------------------------------------------------

// p-normalize collection that is also abs_sym.
// "abs_sym" means the collections's element-wise absolute value is symmetric.

template
<
    typename index_type,
    typename value_type,
    typename vals_id_func_max_collection_type
>
bool p_normalize_vectors_and_trans_abs_sym(

// algorithmic options:
    typename precision_traits<value_type>::scalar p,
    typename precision_traits<value_type>::scalar tolerance,
    std::size_t max_iters,

// input/output data:
    vals_id_func_max_collection_type& vecs,

// output data:
    typename precision_traits<value_type>::scalar* vec_scale,
    std::size_t* iters)
{
    typedef typename precision_traits<value_type>::scalar scalar_type;

    if(
        vecs.num_vecs() != vecs.max_size() ||
        p < 1)
    {
        assert(false);
        return false;
    }

    const index_type num_vecs = vecs.num_vecs();

    std::size_t iters_taken = 0;

    if(0 < num_vecs)
    {
        index_type num_vecs_with_zero_norm = 0;

        std::fill(vec_scale, vec_scale + num_vecs, scalar_type(1));

        const scalar_type zero = 0;

        std::vector<scalar_type> vec_norms;

        try
        {
            vec_norms.resize(num_vecs);
        }
        catch(const std::exception& /*exc*/)
        {
            assert(false);
            return false;
        }

        scalar_type* vec_norms_ptr = &vec_norms.front();

        while(iters_taken < max_iters)
        {
            // First compute the norms and trans-norms.

            p_norm_of_vectors
                <index_type, value_type, vals_id_func_max_collection_type>
                    (p, vecs, vec_norms_ptr);

            // We don't expect there will be any (or many) vectors whose norms
            // will be zeros.  But if there are, we want the algorithm to work
            // without trouble.  Of course, the counts and locations of such
            // vectors will not change as iterations proceed, so we do it once
            // only.

            if(iters_taken == 0)
            {
                num_vecs_with_zero_norm = static_cast<index_type>(std::count(
                    vec_norms.begin(), vec_norms.end(), zero));

                if(num_vecs_with_zero_norm == num_vecs)
                {
                    // all zero matrix, so break.
                    break;
                }
            }

            std::transform(
                vec_norms.begin(), vec_norms.end(), vec_norms.begin(),
                std::ptr_fun<scalar_type, scalar_type>(std::sqrt));

            const scalar_type vec_limit = 1;

            if(tolerance >= 0 &&
                vector_utils_max_abs_diff_with_ignore(
                    num_vecs, vec_norms_ptr, vec_limit, zero) <= tolerance)
            {
                break;
            }

            vector_utils_invert_non_zero(num_vecs, vec_norms_ptr);

            scale_vectors_and_trans
                <index_type, scalar_type, vals_id_func_max_collection_type>
                    (vec_norms_ptr, vecs);

            std::transform(vec_norms.begin(), vec_norms.end(), vec_scale, vec_scale,
                std::multiplies<scalar_type>());

            ++iters_taken;
        }
    }

    if(iters)
        *iters = iters_taken;

    return true;
}

// -----------------------------------------------------------------------------

#endif // P_NORMALIZE_VECTORS_AND_TRANS_H
