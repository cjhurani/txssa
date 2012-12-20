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


#ifndef DENSE_VECTORS_H
#define DENSE_VECTORS_H

// -----------------------------------------------------------------------------

#include "dense_vectors/dense_vectors_utils.h"
#include "math/precision_traits.h"
#include "math/complex_types.h"       // std::abs_square
#include "cpp/const_modifications.h"  // remove_const
#include "internal_api_error/internal_api_error.h"

#include <string>
#include <algorithm>                  // std::swap
#include <cstddef>
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Define simple C++ adaptor classes to wrap sparse/dense vectors
// within a common interface so that template polymorphism works.
// -----------------------------------------------------------------------------

// dense_vectors_id_func supplies actual component ids given an id.  This is
// useful in defining a single interface to iterate over a vector and have a
// knowledge of actual component id.  This is not to be used for increments in
// memory.  It is for a "mathematical" use.  For increments in memory, "inc"
// member functions can be used.

template<typename index_type>
class dense_vectors_id_func
{
public:

    dense_vectors_id_func(index_type in_num_ids)
        :
        num_ids(in_num_ids)
    {
    }

    index_type operator()(index_type i) const
    {
        assert(i < num_ids);
        return i;
    }

private:

    index_type num_ids;
};

// -----------------------------------------------------------------------------

// Dense vectors with uniform offsets (between their "begin").

template<typename index_t, typename value_t>
class dense_vectors
{
public:

    typedef index_t index_type;
    typedef value_t value_type;

    // For dense vectors, id_func_type and inv_id_func_type are the same.
    typedef dense_vectors_id_func<index_type> id_func_type;
    typedef dense_vectors_id_func<index_type> inv_id_func_type;

    dense_vectors()
        :
        n_vecs(0),
        each_vec_size(0),
        leading_dim(0),
        values(0),
        self_allocated(false)
    {
    }

    dense_vectors(
        index_type  in_n_vecs,
        index_type  in_each_vec_size,
        index_type  in_leading_dim,
        value_type* in_values)
        :
        n_vecs(in_n_vecs),
        each_vec_size(in_each_vec_size),
        leading_dim(in_leading_dim),
        values(in_values),
        self_allocated(false)
    {
        // in_values can be 0
        assert(0 < in_leading_dim);
        assert(in_each_vec_size <= in_leading_dim);
    }

    bool use_memory(
        index_type  in_n_vecs,
        index_type  in_each_vec_size,
        index_type  in_leading_dim,
        value_type* in_values)
    {
        // in_values can be 0

        bool success =
            0 < in_leading_dim &&
            in_each_vec_size <= in_leading_dim;

        assert(success);

        if(!success)
            internal_api_error_set_last(
                "dense_vectors<index_t, value_t>::use_memory: Error.");

        dense_vectors tmp(
            in_n_vecs,
            in_each_vec_size,
            in_leading_dim,
            in_values);

        swap(tmp);

        return success;
    }

    bool use_memory(dense_vectors& other)
    {
        return use_memory(
            other.num_vecs(),
            other.vec_size(),
            other.leading_dimension(),
            other.vec_values());
    }

    ~dense_vectors()
    {
        deallocate();
    }

    index_type num_vecs() const
    {
        return n_vecs;
    }

    index_type vec_size() const
    {
        return each_vec_size;
    }

    // Compatibility with sparse vectors
    index_type max_size() const
    {
        return each_vec_size;
    }

    index_type leading_dimension() const
    {
        return leading_dim;
    }

    // Compatibility with sparse vectors
    index_type num_vec_entries(index_type i) const
    {
        (void) i;
        assert(i < n_vecs);

        return each_vec_size;
    }

    const value_type* vec_values() const
    {
        return values;
    }

    value_type* vec_values()
    {
        return values;
    }

    const value_type* vec_values_begin(index_type i) const
    {
        assert(i < n_vecs);

        return values + std::size_t(i) * std::size_t(leading_dim);
    }

    value_type* vec_values_begin(index_type i)
    {
        assert(i < n_vecs);

        return values + std::size_t(i) * std::size_t(leading_dim);
    }

    id_func_type id_func(index_type i) const
    {
        assert(i < n_vecs);

        return id_func_type(each_vec_size);
    }

    inv_id_func_type inv_id_func(index_type i) const
    {
        (void) i;
        assert(i < n_vecs);

        return inv_id_func_type(each_vec_size);
    }

    index_type inc(index_type i) const
    {
        (void) i;
        assert(i < n_vecs);

        return index_type(1);
    }

    void swap(dense_vectors<index_type, value_type>& other)
    {
        std::swap(n_vecs, other.n_vecs);
        std::swap(each_vec_size, other.each_vec_size);
        std::swap(leading_dim, other.leading_dim);
        std::swap(values, other.values);
        std::swap(self_allocated, other.self_allocated);
    }

    void deallocate()
    {
        if(self_allocated)
        {
            delete[] values;
        }
    }

    bool allocate(
        index_type in_n_vecs,
        index_type in_vec_size,
        index_type in_leading_dim = 0)
    {
        if(in_leading_dim < in_vec_size && in_leading_dim != 0)
        {
            assert(false);

            internal_api_error_set_last(
                "dense_vectors<index_t, value_t>::allocate:"
                    " Unacceptable input argument(s).");

            return false;
        }

        // C++ Idiom: Create temporary and swap

        dense_vectors<index_type, value_type> tmp;

        tmp.n_vecs        = in_n_vecs;
        tmp.each_vec_size = in_vec_size;

        if(in_leading_dim == 0)
            tmp.leading_dim = in_vec_size;
        else
            tmp.leading_dim = in_leading_dim;

        // Do this before allocation to allow cleanup in case of exception.
        tmp.self_allocated = true;

        try
        {
            tmp.values = new value_type[
                std::size_t(in_n_vecs) * std::size_t(tmp.leading_dim)];

            deallocate();

            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("dense_vectors<index_t, value_t>::allocate: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    typename remove_const<typename precision_traits<value_type>::scalar>::type
        frobenius_norm_squared() const
    {
        typename precision_traits<value_type>::scalar result = 0;

        const value_type* values_it = values;

        for(index_type j = 0; j < n_vecs; ++j)
        {
            for(index_type i = 0; i < each_vec_size; ++i)
            {
                result += std::abs_square(values_it[i]);
            }

            values_it += leading_dim;
        }

        return result;
    }

    bool axpby(
        const dense_vectors& x,
        value_type a,
        value_type b)
    {
        bool success =
            n_vecs == x.n_vecs &&
            each_vec_size == x.each_vec_size &&
            values &&
            x.values;

        assert(success);

        if(success)
        {
            success = dense_vectors_utils_axpby(
                n_vecs,
                each_vec_size,
                x.values,
                x.leading_dim,
                values,
                leading_dim,
                a,
                b);
        }

      assert(success);

        if(!success)
            internal_api_error_set_last(
                "dense_vectors<index_t, value_t>::axpby:"
                " Error.");

        return success;
    }

    bool add(const dense_vectors& x)
    {
        return axpby(x, value_type(1), value_type(1));
    }

    bool fill(const value_type& fill_with)
    {
        bool success = dense_vectors_utils_fill(
            n_vecs,
            each_vec_size,
            values,
            leading_dim,
            fill_with);

      assert(success);

        if(!success)
            internal_api_error_set_last(
                "dense_vectors<index_t, value_t>::(:"
                " Error.");

        return success;
    }

private:

    index_type  n_vecs;
    index_type  each_vec_size;
    index_type  leading_dim;
    value_type* values;
    bool        self_allocated;

    // Not yet.
    dense_vectors(const dense_vectors&);
    dense_vectors& operator=(const dense_vectors&);
};

// -----------------------------------------------------------------------------

#endif // DENSE_VECTORS_H
