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


#ifndef SPARSE_VECTORS_H
#define SPARSE_VECTORS_H

// -----------------------------------------------------------------------------

#include "sparse_vectors_transpose.h" // sparse_vectors_transpose_nnz_add
#include "math/complex_types.h"       // std::abs_square
#include "math/precision_traits.h"
#include "cpp/std_new_features.h"     // std_new_features_is_sorted
#include "cpp/std_utils.h"            // std_utils_not_in_range, std_utils_not_in_range_closed
#include "internal_api_error/internal_api_error.h"

#include <vector>
#include <algorithm>   // std::{copy, swap, find_if, fill}
#include <numeric>     // std::partial_sum
#include <functional>  // std::{less, greater}
#include <stdexcept>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Classes to either wrap preallocated "compressed sparse vectors"
// with and without values or allocate and own the data.  We use the "vectors
// terminology" and not "row/column terminology", so that these classes work for
// a matrix and its transpose and do not depend on how the matrix is stored, for
// example, the row/col-oriented terminology.  They can work whatever the
// terminology is, but avoiding row/column terminology leads to less confusion.
// -----------------------------------------------------------------------------

// sparse_vectors_id_func supplies actual component ids given an id.
// This is useful in defining a single interface to iterate over a vector and
// have a knowledge of actual component id.  This is not to be used for
// increments in memory.  It is for a "mathematical" use.  For increments in
// memory, "inc" member functions can be used.
//
// sparse_vectors_inv_id_func does the opposite.  Given a component
// id, it returns the offset.

template<typename index_type>
class sparse_vectors_id_func
{
public:

    sparse_vectors_id_func(
        const index_type* in_ids,
        index_type in_num_ids)
        :
        ids(in_ids),
        num_ids(in_num_ids)
    {
        assert(in_ids);
    }

    index_type operator()(index_type i) const
    {
        assert(i < num_ids);
        return ids[i];
    }

private:

    const index_type* ids;
    index_type        num_ids;

    sparse_vectors_id_func& operator=(const sparse_vectors_id_func&);
};

template<typename index_type>
class sparse_vectors_inv_id_func
{
public:

    sparse_vectors_inv_id_func(
        const index_type* in_ids,
        index_type in_num_ids,
        index_type in_max_size)
        :
        ids(in_ids),
        num_ids(in_num_ids),
        max_size(in_max_size)
    {
        assert(in_ids);
        assert(std_new_features_is_sorted(in_ids, in_ids + in_num_ids));
    }

    index_type operator()(index_type i_mapped) const
    {
        assert(i_mapped < max_size);

        // Find i_mapped in ids, and return the offset.
        // Assume ids are sorted.

        const index_type* i_off =
            std::lower_bound(ids, ids + num_ids, i_mapped);

        assert(i_off != ids + num_ids);
        assert(*i_off == i_mapped);

        return index_type(i_off - ids);
    }

private:

    const index_type* ids;
    index_type        num_ids;
    index_type        max_size;

    sparse_vectors_inv_id_func& operator=(const sparse_vectors_inv_id_func&);
};

// -----------------------------------------------------------------------------

// Just the ids, no values.  We specify a different template parameter --
// offset_type -- instead of reusing index_type so that the number of vector
// entries in vector is not limited due to sizeof(index_type), which can be
// as small as 1 if the user wishes.  It also leads to code where it is
// easier to see where type conflicts can occur.

template<typename index_t, typename offset_type>
class sparse_vectors_ids
{
public:

    typedef index_t index_type;

    typedef sparse_vectors_id_func<index_type>     id_func_type;
    typedef sparse_vectors_inv_id_func<index_type> inv_id_func_type;

    sparse_vectors_ids()
        :
        n_vecs(0),
        max_vec_size(0),
        offsets(0),
        ids(0),
        self_allocated(false)
    {
        // This could be made into a compile-time assertion.
        assert(sizeof(index_type) <= sizeof(offset_type));
    }

    sparse_vectors_ids(
        index_type   in_n_vecs,
        index_type   in_max_vec_size,
        offset_type* in_offsets,
        index_type*  in_ids)
        :
        n_vecs(in_n_vecs),
        max_vec_size(in_max_vec_size),
        offsets(in_offsets),
        ids(in_ids),
        self_allocated(false)
    {
        assert(in_offsets);
        assert(in_offsets[0] == 0);
        // assert(in_ids); // in_ids can be 0, so don't assert.
        assert(std_new_features_is_sorted(in_offsets, in_offsets + std::size_t(in_n_vecs) + 1));

        if(in_ids)
        {
            assert(std::find_if(
                in_ids,
                in_ids + in_offsets[in_n_vecs],
                std_utils_not_in_range<index_type>(0, in_max_vec_size)) ==
                in_ids + in_offsets[in_n_vecs]);
        }
    }

    virtual ~sparse_vectors_ids()
    {
        deallocate();
    }

    void use_memory(
        index_type   in_n_vecs,
        index_type   in_max_vec_size,
        offset_type* in_offsets,
        index_type*  in_ids)
    {
        sparse_vectors_ids tmp(
            in_n_vecs,
            in_max_vec_size,
            in_offsets,
            in_ids);

        swap(tmp);
    }

    index_type num_vecs() const
    {
        return n_vecs;
    }

    index_type max_size() const
    {
        return max_vec_size;
    }

    index_type num_vec_entries(index_type i) const
    {
        assert(i < n_vecs);
        assert(offsets);
        assert(offsets[i] <= offsets[std::size_t(i) + 1]);
        assert(offsets[std::size_t(i) + 1] <= offsets[i] + offset_type(max_vec_size));

        return index_type(offsets[std::size_t(i) + 1] - offsets[i]);
    }

    const index_type* vec_ids_begin(index_type i) const
    {
        assert(i < n_vecs);
        assert(ids);
        assert(offsets);

        return ids + offsets[i];
    }

    index_type* vec_ids_begin(index_type i)
    {
        assert(i < n_vecs);
        assert(ids);
        assert(offsets);

        return ids + offsets[i];
    }

    index_type inc(index_type i) const
    {
        assert(i < n_vecs);

        return index_type(1);
    }

    id_func_type id_func(index_type i) const
    {
        assert(i < n_vecs);

        return id_func_type(vec_ids_begin(i), num_vec_entries(i));
    }

    inv_id_func_type inv_id_func(index_type i) const
    {
        assert(i < n_vecs);

        return inv_id_func_type(
            vec_ids_begin(i), num_vec_entries(i), max_vec_size);
    }

    offset_type num_entries() const
    {
        assert(ids);
        assert(offsets);

        return offsets[n_vecs];
    }

    const index_type* vec_ids() const
    {
        return ids;
    }

    index_type* vec_ids()
    {
        return ids;
    }

    const offset_type* vec_offsets() const
    {
        return offsets;
    }

    offset_type* vec_offsets()
    {
        return offsets;
    }

    void swap(sparse_vectors_ids& other)
    {
        std::swap(n_vecs, other.n_vecs);
        std::swap(max_vec_size, other.max_vec_size);
        std::swap(offsets, other.offsets);
        std::swap(ids, other.ids);
        std::swap(self_allocated, other.self_allocated);
    }

    void deallocate()
    {
        if(self_allocated)
        {
            delete[] offsets;
            delete[] ids;
            self_allocated = false;
        }
    }

    // [size_per_vec, size_per_vec + n_vecs) should be valid
    bool allocate(
        index_type in_n_vecs,
        index_type in_max_vec_size,
        const index_type* size_per_vec)
    {
        bool success = size_per_vec != 0;

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors_ids<index_t, offset_t>::allocate:"
                    " Unacceptable input argument(s).");

            return false;
        }

        success = std::find_if(
            size_per_vec,
            size_per_vec + in_n_vecs,
            std_utils_not_in_range_closed<index_type>(0, in_max_vec_size)) ==
            size_per_vec + in_n_vecs;

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors_ids<index_t, offset_t>::allocate:"
                    " some size_per_vec value is not in [0,in_max_vec_size].");

            return false;
        }

        // C++ Idiom: Create temporary and swap

        sparse_vectors_ids tmp;

        tmp.n_vecs       = in_n_vecs;
        tmp.max_vec_size = in_max_vec_size;

        // Do this before allocation to allow cleanup in case of exception.
        tmp.self_allocated = true;

        try
        {
            tmp.offsets = new offset_type[std::size_t(in_n_vecs) + 1];

            tmp.offsets[0] = 0; // size is always at least 1 so [0] element is there.

            // This can lead to integer overflow more easily than explicit loop below.
            //std::partial_sum(size_per_vec, size_per_vec + in_n_vecs, tmp.offsets + 1);

            for(index_type i = 0; i < in_n_vecs; ++i)
                tmp.offsets[std::size_t(i) + 1] = tmp.offsets[i] + size_per_vec[i];

            tmp.ids = new index_type[tmp.offsets[in_n_vecs]];

            deallocate();

            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors_ids<index_t, offset_t>::allocate: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    // [offsets, offsets + n_vecs] should be valid
    bool allocate_using_offsets(
        index_type in_n_vecs,
        index_type in_max_vec_size,
        const offset_type* in_offsets)
    {
        bool success =
            in_offsets &&
            in_offsets[0] == 0;

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors_ids<index_t, offset_t>::allocate_using_offsets:"
                    " Unacceptable input argument(s).");

            return false;
        }

        success = std_new_features_is_sorted(in_offsets, in_offsets + std::size_t(n_vecs) + 1);

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors_ids<index_t, offset_t>::allocate_using_offsets:"
                    " offsets are not sorted in ascending order.");

            return false;
        }

        // C++ Idiom: Create temporary and swap

        sparse_vectors_ids tmp;

        tmp.n_vecs       = in_n_vecs;
        tmp.max_vec_size = in_max_vec_size;

        // Do this before allocation to allow cleanup in case of exception.
        tmp.self_allocated = true;

        try
        {
            tmp.offsets = new offset_type[std::size_t(in_n_vecs) + 1];

            std::copy(in_offsets, in_offsets + std::size_t(in_n_vecs) + 1, tmp.offsets);

            tmp.ids = new index_type[tmp.offsets[in_n_vecs]];

            deallocate();

            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors_ids<index_t, offset_t>::allocate_using_offsets: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    bool create_dense(
        index_type in_n_vecs,
        index_type in_each_vec_size)
    {
        // C++ Idiom: Create temporary and swap

        sparse_vectors_ids tmp;

        tmp.n_vecs       = in_n_vecs;
        tmp.max_vec_size = in_each_vec_size;

        // Do this before allocation to allow cleanup in case of exception.
        tmp.self_allocated = true;

        try
        {
            tmp.offsets = new offset_type[std::size_t(in_n_vecs) + 1];
            tmp.ids     = new index_type[
                offset_type(in_n_vecs) * offset_type(in_each_vec_size)];

            tmp.offsets[0] = 0; // size is always at least 1 so [0] element is there.

            index_type* ids_it = tmp.ids;

            for(index_type i = 0; i < in_n_vecs; ++i, ids_it += in_each_vec_size)
            {
                tmp.offsets[std::size_t(i) + 1] = tmp.offsets[i] + in_each_vec_size;

                for(index_type j = 0; j < in_each_vec_size; ++j)
                {
                    ids_it[j] = j;
                }
            }

            deallocate();

            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors_ids<index_t, offset_t>::create_dense: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

protected:

    index_type   n_vecs;
    index_type   max_vec_size;
    offset_type* offsets;
    index_type*  ids;
    bool         self_allocated;

private:

    // Not yet.
    sparse_vectors_ids(const sparse_vectors_ids&);
    sparse_vectors_ids& operator=(const sparse_vectors_ids&);
};

// -----------------------------------------------------------------------------

// Derive from sparse_vectors_ids to also have value_type's.

template
<
    typename index_t,
    typename offset_type,
    typename value_t
>
class sparse_vectors :
    public sparse_vectors_ids<index_t, offset_type>
{
private:

    typedef sparse_vectors_ids<index_t, offset_type> base_type;

public:

    typedef value_t value_type;
    typedef index_t index_type;

    sparse_vectors()
        :
        base_type(),
        values(0)
    {
    }

    sparse_vectors(
        index_type   in_n_vecs,
        index_type   in_max_vec_size,
        offset_type* in_offsets,
        index_type*  in_ids,
        value_type*  in_values)
        :
        base_type(in_n_vecs, in_max_vec_size, in_offsets, in_ids),
        values(in_values)
    {
        assert(in_values);
    }

    ~sparse_vectors()
    {
        if(base_type::self_allocated)
        {
            delete[] values;
        }
    }

    void use_memory(
        index_type   in_n_vecs,
        index_type   in_max_vec_size,
        offset_type* in_offsets,
        index_type*  in_ids,
        value_type*  in_values)
    {
        sparse_vectors tmp(
            in_n_vecs,
            in_max_vec_size,
            in_offsets,
            in_ids,
            in_values);

        swap(tmp);
    }

    const value_type* vec_values_begin(index_type i) const
    {
        assert(i < base_type::n_vecs);
        assert(base_type::offsets);

        return values + base_type::offsets[i];
    }

    value_type* vec_values_begin(index_type i)
    {
        assert(i < base_type::n_vecs);
        assert(base_type::offsets);

        return values + base_type::offsets[i];
    }

    const value_type* vec_values() const
    {
        return values;
    }

    value_type* vec_values()
    {
        return values;
    }

    void swap(sparse_vectors& other)
    {
        base_type::swap(other);
        std::swap(values, other.values);
    }

    void deallocate()
    {
        if(base_type::self_allocated)
        {
            delete[] values;
        }

        base_type::deallocate();
    }

    bool allocate(
        index_type in_n_vecs,
        index_type in_max_vec_size,
        const index_type* size_per_vec)
    {
        // C++ Idiom: Create temporary and swap

        sparse_vectors tmp;

        bool success =
            tmp.base_type::allocate(in_n_vecs, in_max_vec_size, size_per_vec);

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors<index_t, offset_t, value_t>::allocate: Error in base.");

            return false;
        }

        try
        {
            tmp.values = new value_type[tmp.base_type::offsets[in_n_vecs]];
            deallocate();
            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors<index_t, offset_t, value_t>::allocate: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    bool allocate_using_offsets(
        index_type in_n_vecs,
        index_type in_max_vec_size,
        const offset_type* in_offsets)
    {
        // C++ Idiom: Create temporary and swap

        sparse_vectors tmp;

        bool success =
            tmp.base_type::allocate_using_offsets(in_n_vecs, in_max_vec_size, in_offsets);

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors<index_t, offset_t, value_t>::allocate: Error in base.");

            return false;
        }

        try
        {
            tmp.values = new value_type[tmp.base_type::offsets[in_n_vecs]];
            deallocate();
            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors<index_t, offset_t, value_t>::allocate_using_offsets: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    bool create_dense(
        index_type in_n_vecs,
        index_type in_each_vec_size,
        value_type fill_with)
    {
        // C++ Idiom: Create temporary and swap

        sparse_vectors tmp;

        bool success =
            tmp.base_type::create_dense(in_n_vecs, in_each_vec_size);

        if(!success)
        {
            assert(false);

            internal_api_error_set_last(
                "sparse_vectors<index_t, offset_t, value_t>::create_dense:"
                    " Error in base.");

            return false;
        }

        try
        {
            tmp.values = new value_type[tmp.base_type::offsets[in_n_vecs]];
            std::fill(tmp.values, tmp.values + tmp.base_type::offsets[in_n_vecs],
                fill_with);

            deallocate();
            swap(tmp);
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors<index_t, offset_t, value_t>::create_dense: Exception. ") + exc.what()));

            return false;
        }

        return true;
    }

    // template parameters of this class may be const, so use possibly
    // different parameters for incoming trans variable.
    template
    <
        typename index_type_2,
        typename offset_type_2,
        typename value_type_2
    >
    bool get_transpose(
        sparse_vectors<index_type_2, offset_type_2, value_type_2>& trans) const
    {
        std::vector<index_type_2> transpose_nnz;

        // C++ Idiom: Create temporary and swap
        sparse_vectors<index_type_2, offset_type_2, value_type_2> tmp_trans;

        try
        {
            transpose_nnz.resize(base_type::max_vec_size);

            sparse_vectors_transpose_nnz_add(
                base_type::n_vecs, *this, transpose_nnz.begin());

            // flipped sizes
            bool success = tmp_trans.allocate(base_type::max_vec_size, base_type::n_vecs,
                &transpose_nnz.front());

            if(!success)
            {
                assert(false);

                internal_api_error_set_last(
                    "sparse_vectors<index_t, offset_t, value_t>::get_transpose:"
                        " Error in allocation.");

                return false;
            }
        }
        catch(const std::exception& exc)
        {
            assert(false);

            internal_api_error_set_last(
                (std::string("sparse_vectors<index_t, offset_t, value_t>::get_transpose: Exception. ") + exc.what()));

            return false;
        }

        // Reuse space.  "Rename" first.
        std::vector<index_type_2>& ids_done = transpose_nnz;
        std::fill(ids_done.begin(), ids_done.end(), index_type(0));

        for(index_type_2 i = 0; i < base_type::n_vecs; ++i)
        {
            const index_type num_vec_entries = base_type::num_vec_entries(i);
            const value_type* values         = vec_values_begin(i);
            const index_type* vec            = base_type::vec_ids_begin(i);

            for(index_type_2 j = 0; j < num_vec_entries; ++j)
            {
                const index_type id = vec[j];

                assert(id < base_type::max_vec_size);

                tmp_trans.vec_ids_begin   (id)[ids_done[id]] = i;
                tmp_trans.vec_values_begin(id)[ids_done[id]] = values[j];

                ++ids_done[id];
            }
        }

        trans.swap(tmp_trans);

        return true;
    }

    typename precision_traits<value_type>::scalar
        frobenius_norm_squared() const
    {
        typename precision_traits<value_type>::scalar result = 0;

        for(index_type i = 0; i < base_type::n_vecs; ++i)
        {
            for(offset_type jj = base_type::offsets[i]; jj < base_type::offsets[i+1]; ++jj)
            {
                result += std::abs_square(values[jj]);
            }
        }

        return result;
    }

protected:

    value_type* values;

private:
    // Not yet.
    sparse_vectors(const sparse_vectors&);
    sparse_vectors& operator=(const sparse_vectors&);
};

// -----------------------------------------------------------------------------

#endif // SPARSE_VECTORS_H
