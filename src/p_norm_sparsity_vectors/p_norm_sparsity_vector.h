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


#ifndef P_NORM_SPARSITY_VECTOR_H
#define P_NORM_SPARSITY_VECTOR_H

// -----------------------------------------------------------------------------

#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <sstream>
#include <limits>
#include <algorithm>   // std::{upper_bound, max, sort, copy, fill}
#include <cmath>       // std::{fabs, pow}
#include <cstddef>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute the p-norm based sparsity pattern of a sequence of
// floating point real or complex numbers.
//
// The sequence is specified by a pointer, number of elements, and an increment.
// The output pattern depends on a given ratio in [0,1], "p" value for p-"norm",
// and min_num_nnz, the minimum number of entries to keep in the pattern.
//
// First, the definition of p-"norm".
//
// ||v||_p = p-"norm" of a collection of numbers [v_1, v_2, ..., v_n] =
// (sum_{i=1}^n |v_i|^p)^(1/p)  for 1 <= p < inf
//  sum_{i=1}^n |v_i|^p         for 0 <  p < 1    // Note: no 1/p.
//  max_{i} |v_i|               for p = inf
//  number of non-zero |v_i|    for p = 0
//
// It is a "norm" only when 1 <= p <= inf.  Hence, for a general p, it should be
// called a p-"norm".  We simplify it and don't use the quotes below.
//
// The algorithm maximizes the number of discarded entries while constraining
// the p-norm of discarded entries from above.  If the caller asks, we also keep
// a minimum number of non-zeros (assuming there are that many non-zeros).
//
// Large ratio means more entries are preserved and vice versa.  Ignoring
// min_num_nnz for a moment,
//   p-norm of discarded non-zero entries is <= (1 - ratio) * full p-norm. 
//
// Since we try to discard the maximum number of entries, it means if a
// particular entry is discarded, all entries smaller than that must also be
// discarded.
//
// The main reason we define the inequalities in terms of p-norm of discarded
// entries and not in terms of preserved entries has to do with the p=infinity
// case and ratio < 1.  In this case, preserving just the largest entry would
// mean the p-norm of preserved entries is greater than ratio * full p-norm of
// for ALL ratio < 1.  Thus, ratio becomes useless.  To avoid this, we relate
// 1 - ratio and p-norm of discarded entries.
//
// Another reason for using the norm of discarded entries is when p is large
// (but not necessarily infinite) and its interaction with distribution of
// absolute values of entries.  We expect that the sparsification process makes
// more sense in the scenario when there are few large values (in magnitude) and
// the remaining non-zeros are small.  This is the case for matrices coming from
// hp-FEM.  In the other scenario, when there are many large values but few
// values close to zero, sparsification is less natural.  In very simple terms,
// we're talking about the overall convexity/concavity of the curve of sorted
// non-zero values.  So in the first scenario, with many small entries and few
// large entries, moving a large entry from preserved part to discarded part
// will affect the norm of discarded part much more.  Thus, comparing the norm
// of discarded part with total norm has more "discrimination power", specially
// for high p.  For small p, very close to 1 or even less than it, both norms -
// discarded and preserved, are affected nearly equally.  So our definition is
// motivated by p=inf case as well as interaction of high p case with typical
// distribution we expect.
//
// If 0 < min_num_nnz, we always keep at least min_num_nnz entries while keeping
// the p-norm of discarded entries below the threshold.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

// To compute sum of work_val entries in [work_ids[start], work_ids[end]).

template<typename index_type, typename scalar_type>
scalar_type p_norm_sparsity_vector_sum(
    const scalar_type* work_val,
    const index_type* work_ids,
    index_type start,
    index_type end)
{
    scalar_type sum = 0;

    for(index_type k = start; k < end; ++k)
    {
        sum += work_val[work_ids[k]];

        assert(sum != std::numeric_limits<scalar_type>::infinity());
        assert(sum != std::numeric_limits<scalar_type>::quiet_NaN());
    }

    return sum;
}

// -----------------------------------------------------------------------------

// Two functions to compute sums using different p and scale values in
// [work_val + start, work_val + end) and also update the work_val array.

template<typename index_type, typename scalar_type>
scalar_type p_norm_sparsity_vector_sum_p_update(
    scalar_type p,
    scalar_type* work_val,
    const index_type* work_ids,
    index_type start,
    index_type end)
{
    scalar_type sum = 0;

    for(index_type k = start; k < end; ++k)
    {
        const index_type k_id = work_ids[k];

        work_val[k_id] = std::pow(work_val[k_id], p);

        assert(work_val[k_id] != std::numeric_limits<scalar_type>::infinity());

        sum += work_val[k_id];

        assert(sum != std::numeric_limits<scalar_type>::infinity());
        assert(sum != std::numeric_limits<scalar_type>::quiet_NaN());
    }

    return sum;
}

template<typename index_type, typename scalar_type>
scalar_type p_norm_sparsity_vector_sum_p_scale_update(
    scalar_type p,
    scalar_type scale,
    scalar_type* work_val,
    const index_type* work_ids,
    index_type start,
    index_type end)
{
    scalar_type sum = 0;

    for(index_type k = start; k < end; ++k)
    {
        const index_type k_id = work_ids[k];

        work_val[k_id] = std::pow(work_val[k_id] * scale, p);

        assert(work_val[k_id] != std::numeric_limits<scalar_type>::infinity());

        sum += work_val[k_id];

        assert(sum != std::numeric_limits<scalar_type>::infinity());
        assert(sum != std::numeric_limits<scalar_type>::quiet_NaN());
    }

    return sum;
}

// -----------------------------------------------------------------------------

// For the ordering of index "i" < index "j" if work_val[i] < work_val[j]

template<typename index_type, typename scalar_type>
class p_norm_sparsity_vector_id_less
{
public:

    p_norm_sparsity_vector_id_less(
        const scalar_type* in_work_val,
        const index_type& in_num_ids)
        :
        work_val(in_work_val),
        num_ids(in_num_ids)
    {
        assert(in_work_val);
    }

    bool operator()(const index_type& i, const index_type& j) const
    {
        assert(i < num_ids);
        assert(j < num_ids);
        return work_val[i] < work_val[j];
    }

private:

    const scalar_type* work_val;
    const index_type num_ids;

    p_norm_sparsity_vector_id_less& operator=(const p_norm_sparsity_vector_id_less&);
};

// -----------------------------------------------------------------------------

// For the ordering of index "i" < value "j_val" if work_val[i] < j_val

template<typename index_type, typename scalar_type>
class p_norm_sparsity_vector_val_less
{
public:

    p_norm_sparsity_vector_val_less(
        const scalar_type* in_work_val,
        const index_type& in_num_ids)
        :
        work_val(in_work_val),
        num_ids(in_num_ids)
    {
        assert(in_work_val);
    }

    bool operator()(const index_type& i, const scalar_type& j_val) const
    {
        assert(i < num_ids);
        return work_val[i] < j_val;
    }

    // The two operator() below are defined to avoid the warning
    // Microsoft Visual Studio 8\VC\include\xutility(292) : warning C4244:
    // 'argument' : conversion from 'const double' to 'const int', possible
    // loss of data
    // when the code is compiled in debug mode and visual studio STL
    // runs extra checks and tries to ensure that the range on which
    // upper_bound is used is sorted.  Since scalar_type and index_type
    // are meant to be different in general, this leads to conversion
    // issues.
    bool operator()(const scalar_type& j_val, const index_type& i) const
    {
        return !operator()(i, j_val);
    }

    bool operator()(const index_type& i, const index_type& j) const
    {
        return work_val[i] < work_val[j];
    }

private:

    const scalar_type* work_val;
    const index_type num_ids;

    p_norm_sparsity_vector_val_less& operator=(const p_norm_sparsity_vector_val_less&);
};

// -----------------------------------------------------------------------------

template<typename index_type, typename scalar_type>
void p_norm_sparsity_vector_internal(
    scalar_type ratio,
    scalar_type p,
    index_type min_num_nnz,
    index_type n_entries,
    index_type num_non_zero,
    scalar_type max_abs_val,
    index_type* out_num_nnz,
    index_type* out_ids,
    scalar_type* work_val)
{
    if(num_non_zero == 0)
    {
        // We got all zeros.  This is a special case, and we don't consider it
        // as an error even if min_num_nnz > 0 was given.

        *out_num_nnz = 0;
    }
    else if(min_num_nnz == num_non_zero)
    {
        // We got some some non-zero values, and we're asked to keep the same
        // number of non-zeros.  Simple case, no sorting needed, just keep all.

        *out_num_nnz = num_non_zero;
    }
    else if(min_num_nnz < num_non_zero)
    {
        // The typical case.  This requires sorting.  Here min_num_nnz <
        // num_non_zero.  We sort out_ids but the ordering comes from work_val.
        // work_val does not change.  An id is less than other if the
        // corresponding work value is less.  So use a predicate based
        // comparator.

        std::sort(out_ids, out_ids + num_non_zero,
            p_norm_sparsity_vector_id_less<index_type, scalar_type>(
                work_val, n_entries));

        // Now
        // work_val[out_ids[0]]                is smallest and
        // work_val[out_ids[num_non_zero - 1]] is largest.

        if(p == 0)
        {
            // Preserve (roughly) ratio * num_non_zero entries.  Also take care
            // of min_num_nnz if large enough.  Note: 0 < num_non_zero here.
            // Thus,
            // 0 = ratio => 0 = c
            // 0 < ratio < 1 => 0 < c <= num_non_zero
            // 1 = ratio => c = num_non_zero

            scalar_type c = std::max(scalar_type(1), std::ceil(ratio * scalar_type(num_non_zero)));

            // If ratio * num_non_zero is almost an integer, increase c
            // so that discarded part is "<" and not just "<=" something.

            if(std::fabs(c - ratio * scalar_type(num_non_zero)) <
                100*std::numeric_limits<scalar_type>::epsilon())  // MAGIC CONSTANT
            {
                ++c;
            }

            // c is integral in value but not in type.

            *out_num_nnz = std::min(num_non_zero,
                std::max(min_num_nnz, index_type(c)));
        }
        else if(p == std::numeric_limits<scalar_type>::infinity())
        {
            // Discard values < threshold.

            const scalar_type threshold = (1 - ratio) * max_abs_val;

            // work_val[out_ids[0]] is a/the smallest of all

            if(work_val[out_ids[0]] <= threshold)
            {
                // There is definitely something to discard.

                // In code: *out_num_nnz < num_non_zero and p != out_ids for
                // sure.

                // Binary search with our predicate.  out_ids is not sorted by
                // itself, so binary search does not make sense with the default
                // predicate.

                const index_type* p =
                    std::upper_bound
                        <
                            const index_type*,
                            scalar_type,
                            p_norm_sparsity_vector_val_less<index_type, scalar_type>
                        >(
                        out_ids,
                        out_ids + (num_non_zero - min_num_nnz),
                        threshold,
                        p_norm_sparsity_vector_val_less<index_type, scalar_type>(
                            work_val, num_non_zero));

                // "p - out_ids" values in
                //   [work_val, work_val + (num_non_zero - min_num_nnz)]
                // are "<=" than threshold.  Moreover, "p - out_ids" > 0.
                // Discard them.

                *out_num_nnz = index_type(num_non_zero - index_type(p - out_ids));
            }
            else
            {
                // ratio is large enough that nothing can be discarded
                *out_num_nnz = num_non_zero;
            }
        }
        else
        {
            // For 0 < p < inf cases:
            // Although sparsity pattern is defined in terms of p-norm of
            // discarded entries, for 0 < p < inf, we find the to-be-preserved
            // part and its norm first and deduce the norm of the discarded
            // part.  This is done because of the assumption that in a typical
            // scenario, many more entries will be discarded than preserved
            // (otherwise sparsifying is not so useful).  We want to do
            // computation on fewer entries.  If ratio = 1, we don't have
            // to do any computation.

            if(0 < ratio && ratio < 1 - 100*std::numeric_limits<scalar_type>::epsilon()) // MAGIC CONSTANT
            {
                // range = [mid, num_non_zero) has to be preserved in this case.
                // 0 <= mid because of conditions above.
                // This means [ work_val[out_ids[ range ]] have to be preserved.

                const index_type mid = index_type(num_non_zero - std::max(min_num_nnz, index_type(1)));

                scalar_type val_1, val_2;

                if(p == 1)
                {
                    val_1 = p_norm_sparsity_vector_sum(
                        work_val, out_ids, index_type(0), mid);

                    val_2 = p_norm_sparsity_vector_sum(
                        work_val, out_ids, mid, num_non_zero);
                }
                else if(1 < p)
                {
                    const scalar_type scale = 1/max_abs_val;

                    val_1 = p_norm_sparsity_vector_sum_p_scale_update(
                        p, scale, work_val, out_ids, index_type(0), mid);

                    val_2 = p_norm_sparsity_vector_sum_p_scale_update(
                        p, scale, work_val, out_ids, mid, num_non_zero);
                }
                else // p < 1
                {
                    val_1 = p_norm_sparsity_vector_sum_p_update(p, work_val, out_ids,
                        index_type(0), mid);

                    val_2 = p_norm_sparsity_vector_sum_p_update(p, work_val, out_ids,
                        mid, num_non_zero);
                }

                const scalar_type vec_norm_tmp = val_1 + val_2;

                const scalar_type threshold_to_keep = vec_norm_tmp * (
                    1 < p ? 1 - std::pow(1 - ratio, p) : ratio);

                scalar_type cumulative = val_2;
                index_type k = mid;

                while(0 < k && cumulative < threshold_to_keep)
                {
                    cumulative += work_val[out_ids[--k]];
                }

                *out_num_nnz = index_type(num_non_zero - k);
            }
            else
            {
                if(1 - 100*std::numeric_limits<scalar_type>::epsilon() <= ratio) // MAGIC CONSTANT
                {
                    *out_num_nnz = num_non_zero;
                }
                else // ratio == 0
                {
                    *out_num_nnz = min_num_nnz < index_type(1) ? index_type(1) : min_num_nnz;
                }
            }
        }
    }

    assert(*out_num_nnz <= num_non_zero); // If not, internal error.

    // Shift needed ids to the beginning of output memory.

    if(*out_num_nnz != num_non_zero) // If shifting needs to be done
    {
        assert(0 < *out_num_nnz); // If not, internal error.

        // Note that result iterator does not point to an element in the
        // [first,last), so copy is all right.

        std::copy(
            out_ids + (num_non_zero - *out_num_nnz),
            out_ids + num_non_zero,
            out_ids);

#ifdef _DEBUG

        // Do a little extra work in debug version to fill remaining space with
        // invalid ids.  Invalid in the sense that it is larger than any id
        // we encountered in code earlier but still representable as an
        // index_type.  Thus, using "num_non_zero" is perfect.  This is done to
        // avoid unintentional errors in caller.

        std::fill(
            out_ids + *out_num_nnz,
            out_ids + num_non_zero,
            num_non_zero);
#endif
    }
}

// -----------------------------------------------------------------------------

// Sparsity using p-norm for all p >= 0.

// The output ids are not sorted necessarily.  They also do not depend on value
// of inc.  They are ids computed as if inc were 1.

template<typename index_type, typename value_type>
bool p_norm_sparsity_vector(

// algorithmic options:
    typename precision_traits<value_type>::scalar ratio, // [in] must be in [0,1]
    typename precision_traits<value_type>::scalar p,     // [in] must be in [0,inf]
    index_type min_num_nnz,    // [in] minimum number of non-zeros wanted

// input data:
    index_type n_entries,      // [in] number of entries
    const value_type* v_begin, // [in] beginning of input data
    index_type inc,            // [in] increment to get to next value

// output ids:
    index_type* out_num_nnz,   // [out] number of non-zeros in pattern
    index_type* out_ids,       // [out] size n_entries, use first out_num_nnz

// workspace:                  // [work] size n_entries temporary memory
    typename precision_traits<value_type>::scalar* work_val)
{
    bool success =
        0 <= ratio && ratio <= 1 &&
        0 <= p &&
        min_num_nnz <= n_entries &&
        v_begin &&
        0 < inc &&
        out_num_nnz &&
        out_ids &&
        work_val;

    if(!success)
    {
        std::ostringstream oss;
        oss << "p_norm_sparsity_vector: Unacceptable input argument(s).";
        if(!(0 <= ratio && ratio <= 1)) oss << " !(0 <= ratio && ratio <= 1).";
        if(!(0 <= p))                   oss << " !(0 <= p).";
        if(!(min_num_nnz <= n_entries)) oss << " !(min_num_nnz <= n_entries)).";
        if(!(v_begin))                  oss << " !(v_begin).";
        if(!(0 < inc))                  oss << " !(0 < inc).";
        if(!(out_num_nnz))              oss << " !(out_num_nnz).";
        if(!(out_ids))                  oss << " !(out_ids).";
        if(!(work_val))                 oss << " !(work_val).";

        internal_api_error_set_last(oss.str());
        return false;
    }

    typedef typename precision_traits<value_type>::scalar scalar_type;

    scalar_type max_abs_val = -1;
    index_type num_non_zero = 0;

    // Count number of non-zeros and fill workspace and find maximum
    // absolute value.  max_abs_val is needed only when 1 < p.  For
    // 1 < p, we use it for scaling down values before norm computation
    // and it is used as infinity-norm for for p == inf.

    const value_type* v_it = v_begin; 

    if(1 < p)
    {
        max_abs_val = 0;

        for(index_type j = 0; j < n_entries; ++j)
        {
            const value_type val = *v_it;
            v_it += inc;

            if(val != value_type(0))
            {
                const scalar_type abs_val = std::abs(val);

                if(max_abs_val < abs_val)
                    max_abs_val = abs_val;

                work_val[j] = abs_val;
                out_ids[num_non_zero++] = j;
            }
        }
    }
    else
    {
        for(index_type j = 0; j < n_entries; ++j)
        {
            const value_type val = *v_it;
            v_it += inc;

            if(val != value_type(0))
            {
                work_val[j] = std::abs(val);
                out_ids[num_non_zero++] = j;
            }
        }
    }

    if(num_non_zero == 0 || min_num_nnz <= num_non_zero)
    {
        // In the debug version, fill the extra space in out_ids with a
        // large value depending on index_type.

#ifdef _DEBUG

        std::fill(
            out_ids + num_non_zero,
            out_ids + n_entries,
            std::numeric_limits<index_type>::max());
#endif

        // If everything above worked, this function should always work.

        p_norm_sparsity_vector_internal(ratio, p, min_num_nnz, n_entries,
            num_non_zero, max_abs_val, out_num_nnz, out_ids, work_val);
    }
    else
    {
        // 0 < num_non_zero < min_num_nnz.
        success = false;
        assert(false);

        const char* error_string =
            "p_norm_sparsity_vector: Unacceptable input value(s)."
            " 0 < num_non_zero < min_num_nnz.";
        internal_api_error_set_last(error_string);
    }

#if defined(_DEBUG) || defined(P_NORM_SPARSITY_VECTOR_TEST)

    if(success)
    {
        const char* error_string =
            "p_norm_sparsity_vector: Internal error.";

        // Verify that there is no internal error.

        if(*out_num_nnz < min_num_nnz)
        {
            assert(false);
            internal_api_error_set_last(error_string);
            return false;
        }

        scalar_type kept_norm = 0;

        if(p == std::numeric_limits<scalar_type>::infinity())
        {
            for(index_type j = 0; j < *out_num_nnz; ++j)
            {
                const value_type val = v_begin[std::size_t(inc) * std::size_t(out_ids[j])];

                if(val != value_type(0))
                {
                    const scalar_type abs_val = std::abs(val);

                    if(kept_norm < abs_val)
                        kept_norm = abs_val;
                }
                else
                {
                    assert(false);
                    internal_api_error_set_last(error_string);
                    return false;
                }
            }

            if(kept_norm != max_abs_val)
            {
                assert(false);
                internal_api_error_set_last(error_string);
                return false;
            }
        }
        else
        {
            for(index_type j = 0; j < *out_num_nnz; ++j)
            {
                const value_type val = v_begin[std::size_t(inc)*std::size_t(out_ids[j])];

                if(val != value_type(0))
                {
                    kept_norm += std::pow(std::abs(val), p);
                }
                else
                {
                    assert(false);
                    internal_api_error_set_last(error_string);
                    return false;
                }
            }

            const value_type* v_it = v_begin; 

            scalar_type full_norm = 0;

            for(index_type j = 0; j < n_entries; ++j)
            {
                const value_type val = *v_it;
                v_it += inc;

                if(val != value_type(0))
                {
                    full_norm += std::pow(std::abs(val), p);
                }
            }

            scalar_type discarded_norm;
            scalar_type threshold;

            if(p <= 1)
            {
                discarded_norm = full_norm - kept_norm;
                threshold = (1 - ratio)*full_norm;
            }
            else
            {
                if(full_norm < kept_norm)
                {
                    if(kept_norm - full_norm > 100*std::numeric_limits<scalar_type>::epsilon()*full_norm) // MAGIC CONSTANT
                    {
                        assert(false);
                        internal_api_error_set_last(error_string);
                        return false;
                    }
                }

                discarded_norm = std::pow(full_norm - kept_norm, 1/p);
                threshold = (1 - ratio)*std::pow(full_norm, 1/p);
            }

            if(ratio == 1 && *out_num_nnz != num_non_zero)
            {
                assert(false);
                internal_api_error_set_last(error_string);
                return false;
            }

            if(ratio < 1 && discarded_norm >= threshold)
            {
                assert(false);
                internal_api_error_set_last(error_string);
                return false;
            }
        }
    }

#endif

    return success;
}

// -----------------------------------------------------------------------------

#endif // P_NORM_SPARSITY_VECTOR_H
