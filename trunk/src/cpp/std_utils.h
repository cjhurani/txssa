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


#ifndef STD_UTILS_H
#define STD_UTILS_H

// -----------------------------------------------------------------------------

#include <complex>     // std::abs for complex
#include <cstddef>
#include <cmath>       // std::abs for real, std::fabs
#include <limits>
#include <iterator>
#include <algorithm>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Utility functions/classes similar to C++ standard library.
// -----------------------------------------------------------------------------

template<typename value_type>
class std_utils_not_in_range
{
public:

    // Range is [in_v1, in_v2)

    std_utils_not_in_range(const value_type& in_v1, const value_type& in_v2)
        :
        v1(in_v1),
        v2(in_v2)
    {
        assert(in_v1 < in_v2);
    }

    bool operator()(const value_type& v) const
    {
        // Use !(v < v2) instead of v2 <= v2, so that only "<" is used.

        return v < v1 || !(v < v2);
    }

private:

    const value_type& v1;
    const value_type& v2;

    std_utils_not_in_range& operator=(const std_utils_not_in_range&);
};

// -----------------------------------------------------------------------------

template<typename value_type>
class std_utils_not_in_range_closed
{
public:

    // Range is [in_v1, in_v2]

    std_utils_not_in_range_closed(const value_type& in_v1, const value_type& in_v2)
        :
        v1(in_v1),
        v2(in_v2)
    {
        assert(in_v1 < in_v2);
    }

    bool operator()(const value_type& v) const
    {
        return v < v1 || v2 < v;
    }

private:

    const value_type& v1;
    const value_type& v2;

    std_utils_not_in_range_closed& operator=(const std_utils_not_in_range_closed&);
};

// -----------------------------------------------------------------------------

// Count of the number of values whose abs <= threshold starting from the end.
// Stop counting if any value is encountered whose abs is greater than
// threshold.

template
<
    typename index_type,
    typename value_type,
    typename abs_value_type
>
index_type std_utils_count_less_equal_abs_reverse_inc(
    const value_type* begin,
    index_type num_vals,
    std::size_t inc,
    abs_value_type threshold)
{
    assert(0 < inc);
    assert(begin);

    index_type count = 0;

    // avoid over-flow by using size_t.
    const value_type* end = begin + std::size_t(num_vals) * inc;

    while(end != begin)
    {
        end -= inc;

        if(std::abs(*end) <= threshold)
            ++count;
        else
            break;
    }

    return count;
}

// -----------------------------------------------------------------------------

// The objective is to compute minimum and maximum values in a sequence but when
// the sequence is seen as divided into two at "separated_at".  Thus, two max
// and two min values are found.  The separation leads to a closed half interval
// junction for both left and right.  Only "<" and "==" comparison operators are
// used for value_t.

template<typename value_t>
struct std_utils_separated_min_max
{
    typedef value_t value_type;

    std_utils_separated_min_max(value_type in_separated_at = value_type())
        :
        max_r(-std::numeric_limits<value_type>::max()),
        min_r( std::numeric_limits<value_type>::max()),
        max_l(-std::numeric_limits<value_type>::max()),
        min_l( std::numeric_limits<value_type>::max()),
        separated_at(in_separated_at)
    {
    }

    void use(value_type v)
    {
        if(v == separated_at)
        {
            max_l = separated_at;
            min_r = separated_at;
        }
        else if(v < separated_at)
        {
            if(max_l < v)
                max_l = v;

            if(v < min_l)
                min_l = v;
        }
        else // if(separated_at < v)
        {
            if(max_r < v)
                max_r = v;

            if(v < min_r)
                min_r = v;
        }
    }

    // Was any value seen on the left side (including separated_at)
    bool any_in_left() const
    {
        return !(max_l < min_l);
    }

    // Was any value seen on the right side (including separated_at)
    bool any_in_right() const
    {
        return !(max_r < min_r);
    }

    // Was any value seen on the left side (not including separated_at)
    bool any_in_strict_left() const
    {
        return min_l < separated_at;
    }

    // Was any value seen on the right side (not including separated_at)
    bool any_in_strict_right() const
    {
        return separated_at < max_r;
    }

    bool any_at_separation() const
    {
        return max_l == min_r;
    }

    const value_type& max_left() const
    {
        return max_l;
    }

    const value_type& min_left() const
    {
        return min_l;
    }

    const value_type& max_right() const
    {
        return max_r;
    }

    const value_type& min_right() const
    {
        return min_r;
    }

    const value_type& separation() const
    {
        return separated_at;
    }

    // Return true if did perturb.
    bool perturb(value_type fuzz = value_type(1E2)) // MAGIC CONSTANT
    {
        bool perturbed = false;

        // If the values are to be perturbed
        if(any_in_strict_left() && any_in_strict_right())
        {
            const value_type max_left_dist  = separated_at - max_l;
            const value_type min_left_dist  = separated_at - min_l;
            const value_type max_right_dist = max_r - separated_at;
            const value_type min_right_dist = min_r - separated_at;

            const value_type large_diff = std::fabs(max_right_dist - min_left_dist);
            const value_type small_diff = std::fabs(min_right_dist - max_left_dist);

            const value_type rel_tol = fuzz * std::numeric_limits<value_type>::epsilon();

            const value_type large_tol = rel_tol * (max_right_dist + min_left_dist);
            const value_type small_tol = rel_tol * (min_right_dist + max_left_dist);

            if(large_diff <= large_tol && small_diff <= small_tol)
            {
                // We have a sequence such that when viewed on both sides of
                // separated_at, the values are roughly within same min/max
                // range.  This does not mean that for every left value there
                // is one and only one right value.  But it makes sense that we
                // should produce a binning that tries to adjust somewhat for
                // this scenario even in presence of floating point errors in
                // incoming data.

                // Choose value that is closer to separated_at between max_left, min_right
                if(max_left_dist < min_right_dist)
                    min_r = separated_at + max_left_dist;
                else
                    max_l = separated_at - min_right_dist;

                // Choose value that is farther from separated_at between max_right, min_left
                if(min_left_dist < max_right_dist)
                    min_l = separated_at - max_right_dist;
                else
                    max_r = separated_at + min_left_dist;

                // If the code above works, these asserts should not fail even in
                // presence of any floating point inaccuracies.
                assert((max_r - separated_at) == (separated_at - min_l));
                assert((min_r - separated_at) == (separated_at - max_l));

                perturbed = true;
            }
        }

        return perturbed;
    }

private:

    value_type max_r;
    value_type min_r;
    value_type max_l;
    value_type min_l;
    value_type separated_at;
};

// -----------------------------------------------------------------------------

// Given a range [ids_begin, ids_end) specified by forward iterators, it
// one-to-one maps the values so that the new range of values is smallest
// possible.  In other words, if there are n unique values present in
// [ids_begin, ids_end), then after this mapping, they are all mapped to [0, n).
//
// We don't specifically check for the following conditions but they must hold.
//
// ids_it_fwd must be a forward iterator type.
// work_it_rand must be a random access iterator type.
//
// ids_it_fwd's value_type must be an integral type.
// work_it_rand's value_type must be same as ids_it_fwd's value_type.
// all values in [ids_begin, ids_end) must be less than max_val.
// [work_begin, work_begin + max_val) must be a valid range (and
// it will be over-written).

template<typename ids_it_fwd, typename work_it_rand>
void std_utils_bin_mapping(

// input data:
    typename std::iterator_traits<ids_it_fwd>::value_type max_val,

// input/output data:
    ids_it_fwd ids_begin,
    ids_it_fwd ids_end,

// output:
    typename std::iterator_traits<ids_it_fwd>::value_type& n_unique_vals,

// workspace:
    work_it_rand work_begin)
{
    typedef typename std::iterator_traits<ids_it_fwd>::value_type index_type;

    // Since all ids in [ids_begin, ids_end) are less than max_val, max_val is
    // an impossible_id.
    const index_type impossible_id = max_val;

    std::fill(work_begin, work_begin + max_val, impossible_id);

    n_unique_vals = index_type();

    ids_it_fwd ids_it = ids_begin;

    // Count number of unique values and fill work space with the mapped ids.

    while(ids_it != ids_end)
    {
        assert(*ids_it < max_val);

        index_type& id = work_begin[*ids_it++];

        if(id == impossible_id)
        {
            id = n_unique_vals++;
        }
    }

    ids_it = ids_begin;

    // Copy mapped ids from work space to [ids_begin, ids_end).

    while(ids_it != ids_end)
    {
        *ids_it = work_begin[*ids_it];
        ++ids_it;
    }
}

// -----------------------------------------------------------------------------

#endif // STD_UTILS_H
