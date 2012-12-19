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


#ifndef P_NORM_OF_VECTORS_H
#define P_NORM_OF_VECTORS_H

// -----------------------------------------------------------------------------

#include "math/precision_traits.h"
#include "math/vector_utils.h"
#include "cpp/const_modifications.h"
#include <cmath>       // std::{pow, fabs}
#include <limits>
#include <algorithm>   // std::fill
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute p-norms of each vector in a collection of vectors.  The
// interface may look strange but that is because it is taken from
// p_norm_of_vectors_and_trans.h and modified.
// -----------------------------------------------------------------------------

// Template-concept: vals_collection_type
// Following member functions should be present.
//
// index_type          num_vecs        ()           const
// [const] value_type* vec_values_begin(index_type) const
// index_type          num_vec_entries (index_type) const

// -----------------------------------------------------------------------------

// p_case_value:
// p = 0       =>  0
// 0 < p < 1   =>  -1
// p = 1       =>  1
// 1 < p < inf =>  2
// p = inf     =>  3

template
<
    typename index_type,
    typename value_type,
    int p_case_value
>
class p_norm_of_vectors_helper;

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
class p_norm_of_vectors_helper<index_type, value_type, -1>
{
public:

    typedef typename precision_traits<value_type>::scalar scalar_type;

    p_norm_of_vectors_helper(
        scalar_type* in_vec_norms,
        scalar_type  in_p)
        :
        vec_norms(in_vec_norms),
        p(in_p)
    {
    }

    void process(index_type i, const value_type& v)
    {
        vec_norms[i] += std::pow(std::abs(v), p);
    }

private:

    scalar_type* vec_norms;
    scalar_type  p;
};

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
class p_norm_of_vectors_helper<index_type, value_type, 0>
{
public:

    typedef typename precision_traits<value_type>::scalar scalar_type;

    p_norm_of_vectors_helper(
        scalar_type* in_vec_norms)
        :
        vec_norms(in_vec_norms)
    {
    }

    void process(index_type i, const value_type& v)
    {
        if(v != value_type())
        {
            ++vec_norms[i];
        }
    }

private:

    scalar_type* vec_norms;
};

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
class p_norm_of_vectors_helper<index_type, value_type, 1>
{
public:

    typedef typename precision_traits<value_type>::scalar scalar_type;

    p_norm_of_vectors_helper(
        scalar_type* in_vec_norms)
        :
        vec_norms(in_vec_norms)
    {
    }

    void process(index_type i, const value_type& v)
    {
        vec_norms[i] += std::abs(v);
    }

private:

    scalar_type* vec_norms;
};

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
class p_norm_of_vectors_helper<index_type, value_type, 2>
{
public:

    typedef typename precision_traits<value_type>::scalar scalar_type;

    p_norm_of_vectors_helper(
        scalar_type* in_vec_norms,
        scalar_type  in_p)
        :
        vec_norms(in_vec_norms),
        p(in_p)
    {
    }

    void process(index_type i, const value_type& v)
    {
        vec_norms[i] += std::pow(std::abs(v), p);
    }

private:

    scalar_type* vec_norms;
    scalar_type  p;
};

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
class p_norm_of_vectors_helper<index_type, value_type, 3>
{
public:

    typedef typename precision_traits<value_type>::scalar scalar_type;

    p_norm_of_vectors_helper(
        scalar_type* in_vec_norms)
        :
        vec_norms(in_vec_norms)
    {
    }

    void process(index_type i, const value_type& v)
    {
        const scalar_type to_use = std::abs(v);

        if( vec_norms[i] < to_use)
            vec_norms[i] = to_use;
    }

private:

    scalar_type* vec_norms;
};

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename value_type,
    typename vals_collection_type,
    typename helper_type
>
void p_norm_of_vectors_accumulate(
    const vals_collection_type& vecs,
    helper_type& helper)
{
    const index_type num_vecs = vecs.num_vecs();

    for(index_type i = 0; i < num_vecs; ++i)
    {
        const index_type num_vec_entries = vecs.num_vec_entries(i);

        const value_type* vals = vecs.vec_values_begin(i);

        for(index_type j = 0; j < num_vec_entries; ++j)
        {
            helper.process(i, vals[j]);
        }
    }
}

// -----------------------------------------------------------------------------

template
<
    typename scalar_type,
    typename vals_collection_type
>
void p_norm_of_vectors(
    scalar_type p,
    const vals_collection_type& vecs,
    scalar_type* func_vec_norms)
{
    typedef typename vals_collection_type::index_type index_type1;
    typedef typename remove_const<index_type1>::type index_type;
    typedef typename vals_collection_type::value_type value_type;

    const index_type num_vecs = vecs.num_vecs();

    assert(0 <= p);
    assert(func_vec_norms);
    assert(num_vecs == vecs.num_vecs());

    std::fill(func_vec_norms, func_vec_norms + num_vecs, scalar_type(0));

    const scalar_type inf = std::numeric_limits<scalar_type>::infinity();

    if(p == 0)
    {
        typedef p_norm_of_vectors_helper<index_type, value_type, 0>
            helper_type;

        helper_type helper(func_vec_norms);

        p_norm_of_vectors_accumulate
            <index_type, value_type, vals_collection_type, helper_type>
                (vecs, helper);
    }
    else if(p == 1)
    {
        typedef p_norm_of_vectors_helper<index_type, value_type, 1>
            helper_type;

        helper_type helper(func_vec_norms);

        p_norm_of_vectors_accumulate
            <index_type, value_type, vals_collection_type, helper_type>
                (vecs, helper);
    }
    else if(p == inf)
    {
        typedef p_norm_of_vectors_helper<index_type, value_type, 3>
            helper_type;

        helper_type helper(func_vec_norms);

        p_norm_of_vectors_accumulate
            <index_type, value_type, vals_collection_type, helper_type>
                (vecs, helper);
    }
    else if(p < 1)
    {
        typedef p_norm_of_vectors_helper<index_type, value_type, -1>
            helper_type;

        helper_type helper(func_vec_norms, p);

        p_norm_of_vectors_accumulate
            <index_type, value_type, vals_collection_type, helper_type>
                (vecs, helper);
    }
    else // if(1 < p && p < inf)
    {
        typedef p_norm_of_vectors_helper<index_type, value_type, 2>
            helper_type;

        helper_type helper(func_vec_norms, p);

        p_norm_of_vectors_accumulate
            <index_type, value_type, vals_collection_type, helper_type>
                (vecs, helper);
    }

    // Compute the correct norm from the accumulated values.

    if(1 < p && p < inf)
    {
        vector_utils_replace_with_pow(num_vecs, func_vec_norms, 1/p);
    }
}

// -----------------------------------------------------------------------------

#endif // P_NORM_OF_VECTORS_H
