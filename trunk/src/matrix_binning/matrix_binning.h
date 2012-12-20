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


#ifndef MATRIX_BINNING_H
#define MATRIX_BINNING_H

// -----------------------------------------------------------------------------

#include "cpp/std_utils.h"
#include "cpp/const_modifications.h"
#include "math/precision_traits.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Create a binned pattern depending on matrix values.
// -----------------------------------------------------------------------------

template<typename bin_index_type, typename value_type>
struct matrix_binning_worker
{
    matrix_binning_worker()
        :
        max_n_left_bins(0),
        max_n_right_bins(0),
        inv_h_l(0),
        inv_h_r(0),
        separated_at(0),
        min_left(0),
        max_right(0),
        left_tol(-1),
        right_tol(-1)
    {
    }

    void use(
        const std_utils_separated_min_max<value_type>& sep_min_max,
        bin_index_type max_num_bins)
    {
        // Since we do have some values (otherwise this function would not be
        // called), strictly speaking we have 7 distinct cases of whether we
        // have (or don't have) values in left, right, or separated_at
        // (2^3 - 1 = 7).  We want to bin in a symmetric form of half-open
        // intervals.  Like this
        //
        //    [) [) [) | (] (] (] (] (]
        //    __left__ | _____right____
        //             |
        //        separated_at
        //
        // We want to do this for any kind of input (symmetric or not,
        // present in both sides or not and so on).  The separated_at
        // value will go in its own special bin.

        const bool any_at_separation = sep_min_max.any_at_separation();
        separated_at = sep_min_max.separation();

        const bin_index_type loc_max_num_left_right_bins =
            max_num_bins - bin_index_type(any_at_separation ? 1 : 0);

        min_left  = sep_min_max.min_left();
        max_right = sep_min_max.max_right();

        const value_type left_dist  = separated_at > min_left ? separated_at - min_left : 0;
        const value_type right_dist = max_right > separated_at ? max_right - separated_at : 0;
        const value_type total_dist = left_dist + right_dist;

        max_n_left_bins  = bin_index_type((value_type(loc_max_num_left_right_bins) * left_dist ) / total_dist);
        max_n_right_bins = bin_index_type((value_type(loc_max_num_left_right_bins) * right_dist) / total_dist);

        assert(max_n_left_bins + max_n_right_bins <= loc_max_num_left_right_bins);

        inv_h_l = value_type(max_n_left_bins)/left_dist;
        inv_h_r = value_type(max_n_right_bins)/right_dist;

        left_tol =  100 * std::numeric_limits<value_type>::epsilon() * left_dist;   // MAGIC CONSTANT
        right_tol = 100 * std::numeric_limits<value_type>::epsilon() * right_dist;  // MAGIC CONSTANT
    }

    bin_index_type bin_of(const value_type& v) const
    {
        assert(!(v < min_left));
        assert(!(max_right < v));
        assert(!(left_tol < 0));
        assert(!(right_tol < 0));

        bin_index_type ans;

        if(v == separated_at)
        {
            ans = max_n_left_bins + max_n_right_bins;
        }
        else if(v < separated_at)
        {
            if(separated_at - v <= left_tol)
                ans = max_n_left_bins + max_n_right_bins;
            else
                ans = bin_index_type((v - min_left) * inv_h_l);
        }
        else // if(separated_at < v)
        {
            if(v - separated_at <= right_tol)
                ans = max_n_left_bins + max_n_right_bins;
            else
                ans = bin_index_type((max_right - v) * inv_h_r) + max_n_left_bins;
        }

        return ans;
    }

private:

    bin_index_type max_n_left_bins;
    bin_index_type max_n_right_bins;

    value_type inv_h_l;
    value_type inv_h_r;

    value_type separated_at;
    value_type min_left;
    value_type max_right;

    value_type left_tol;
    value_type right_tol;
};

// -----------------------------------------------------------------------------

template
<
    typename bin_index_type,
    typename ids_collection_type,
    typename vals_inc_collection_type
>
bool matrix_binning(
    const ids_collection_type& sparse_pat,
    const vals_inc_collection_type& vecs,
    bin_index_type max_num_bins,
    bin_index_type* bin_ids,    // size = number of entries in sparse_pat.
    bin_index_type& actual_num_bins,
    bin_index_type* work_array) // size = max_num_bins, used only if max_num_bins > 1
{
    typedef typename vals_inc_collection_type::index_type index_type;

    const index_type num_vecs = vecs.num_vecs();

    bool success =
        bin_ids &&
        (num_vecs == sparse_pat.num_vecs()) &&
        (work_array || max_num_bins <= 1);

    assert(success);

    if(success)
    {
        success = false;

        if(max_num_bins == 0)
        {
            actual_num_bins = 0;

            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *bin_ids++ = actual_num_bins++;

            success = true;
        }
        else if(max_num_bins == 1)
        {
            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *bin_ids++ = bin_index_type();

            actual_num_bins = 1;
            success = true;
        }
        else if(bin_index_type(1) < max_num_bins)
        {
            typedef typename vals_inc_collection_type::value_type value_type_1;
            typedef typename remove_const<value_type_1>::type value_type;
            typedef typename precision_traits<value_type>::scalar scalar_type;

            const scalar_type separated_at = 0;

            std_utils_separated_min_max<value_type> sep_min_max(separated_at);

            for(index_type i = 0; i < num_vecs; ++i)
            {
                const index_type num_sparse_entries = sparse_pat.num_vec_entries(i);
                const index_type* sparsity_ids = sparse_pat.vec_ids_begin(i);
                const value_type* vals = vecs.vec_values_begin(i);
                const index_type inc = vecs.inc(i);

                assert(sparsity_ids);
                assert(num_sparse_entries <= vecs.num_vec_entries(i));
                assert(vals);

                for(index_type j = 0; j < num_sparse_entries; ++j)
                {
                    const value_type v = vals[inc * sparsity_ids[j]];

                    sep_min_max.use(v);
                }
            }

            actual_num_bins = 0;

            // If something to do
            if(sep_min_max.any_in_left() || sep_min_max.any_in_right())
            {
                sep_min_max.perturb();

                matrix_binning_worker<bin_index_type, value_type> worker;
                worker.use(sep_min_max, max_num_bins);

                bin_index_type* bin_ids_it = bin_ids;

                // Go through the values again.
                for(index_type i = 0; i < num_vecs; ++i)
                {
                    const index_type num_sparse_entries = sparse_pat.num_vec_entries(i);
                    const index_type* sparsity_ids = sparse_pat.vec_ids_begin(i);
                    const value_type* vals = vecs.vec_values_begin(i);
                    const index_type inc = vecs.inc(i);

                    for(index_type j = 0; j < num_sparse_entries; ++j)
                    {
                        *bin_ids_it++ = worker.bin_of(vals[inc * sparsity_ids[j]]);
                    }
                }

                std_utils_bin_mapping(
                    max_num_bins,
                    bin_ids, bin_ids_it,
                    actual_num_bins,
                    work_array);
            }

            success = true;
        }
    }

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "matrix_binning: Error.");

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename bin_index_type,
    typename ids_collection_type,
    typename vals_inc_collection_type
>
bool matrix_binning(
    const ids_collection_type& sparse_pat,
    const vals_inc_collection_type& vecs,
    bin_index_type max_num_bins,  // individual max in real and imag.
    bin_index_type* real_bin_ids, // size = number of entries in sparse_pat.
    bin_index_type* imag_bin_ids, // size = number of entries in sparse_pat.
    bin_index_type& real_actual_num_bins,
    bin_index_type& imag_actual_num_bins,
    bin_index_type* work_array) // size = max_num_bins, used only if max_num_bins > 1
{
    typedef typename vals_inc_collection_type::index_type index_type;

    const index_type num_vecs = vecs.num_vecs();

    bool success =
        real_bin_ids &&
        imag_bin_ids &&
        (num_vecs == sparse_pat.num_vecs()) &&
        (work_array || max_num_bins <= 1);

    assert(success);

    if(success)
    {
        success = false;

        if(max_num_bins == 0)
        {
            real_actual_num_bins = 0;
            imag_actual_num_bins = 0;

            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *real_bin_ids++ = real_actual_num_bins++;

            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *imag_bin_ids++ = imag_actual_num_bins++;

            success = true;
        }
        else if(max_num_bins == 1)
        {
            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *real_bin_ids++ = bin_index_type();

            for(index_type i = 0; i < num_vecs; ++i)
                for(index_type j = 0; j < sparse_pat.num_vec_entries(i); ++j)
                    *imag_bin_ids++ = bin_index_type();

            real_actual_num_bins = 1;
            imag_actual_num_bins = 1;
            success = true;
        }
        else if(bin_index_type(1) < max_num_bins)
        {
            typedef typename vals_inc_collection_type::value_type value_type_1;
            typedef typename remove_const<value_type_1>::type value_type;
            typedef typename precision_traits<value_type>::scalar scalar_type;

            const scalar_type separated_at = 0;

            std_utils_separated_min_max<scalar_type> real_sep_min_max(separated_at);
            std_utils_separated_min_max<scalar_type> imag_sep_min_max(separated_at);

            for(index_type i = 0; i < num_vecs; ++i)
            {
                const index_type num_sparse_entries = sparse_pat.num_vec_entries(i);
                const index_type* sparsity_ids = sparse_pat.vec_ids_begin(i);
                const value_type* vals = vecs.vec_values_begin(i);
                const index_type inc = vecs.inc(i);

                assert(sparsity_ids);
                assert(num_sparse_entries <= vecs.num_vec_entries(i));
                assert(vals);

                for(index_type j = 0; j < num_sparse_entries; ++j)
                {
                    const value_type v = vals[inc * sparsity_ids[j]];

                    real_sep_min_max.use(v.real());
                    imag_sep_min_max.use(v.imag());
                }
            }

            real_actual_num_bins = 0;
            imag_actual_num_bins = 0;

            // If something to do
            if(real_sep_min_max.any_in_left() || real_sep_min_max.any_in_right() ||
               imag_sep_min_max.any_in_left() || imag_sep_min_max.any_in_right())
            {
                real_sep_min_max.perturb();
                imag_sep_min_max.perturb();

                matrix_binning_worker<bin_index_type, scalar_type> rworker, iworker;
                rworker.use(real_sep_min_max, max_num_bins);
                iworker.use(imag_sep_min_max, max_num_bins);

                bin_index_type* real_bin_ids_it = real_bin_ids;
                bin_index_type* imag_bin_ids_it = imag_bin_ids;

                // Go through the values again.
                for(index_type i = 0; i < num_vecs; ++i)
                {
                    const index_type num_sparse_entries = sparse_pat.num_vec_entries(i);
                    const index_type* sparsity_ids = sparse_pat.vec_ids_begin(i);
                    const value_type* vals = vecs.vec_values_begin(i);
                    const index_type inc = vecs.inc(i);

                    for(index_type j = 0; j < num_sparse_entries; ++j)
                    {
                        const value_type v = vals[inc * sparsity_ids[j]];

                        *real_bin_ids_it++ = rworker.bin_of(v.real());
                        *imag_bin_ids_it++ = iworker.bin_of(v.imag());
                    }
                }

                bin_index_type* real_bin_ids_end = real_bin_ids_it;
                bin_index_type* imag_bin_ids_end = imag_bin_ids_it;

                std_utils_bin_mapping(
                    max_num_bins,
                    real_bin_ids, real_bin_ids_end,
                    real_actual_num_bins,
                    work_array);

                std_utils_bin_mapping(
                    max_num_bins,
                    imag_bin_ids, imag_bin_ids_end,
                    imag_actual_num_bins,
                    work_array);
            }

            success = true;
        }
    }

    assert(success);

    if(!success)
        internal_api_error_set_last(
            "matrix_binning: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // MATRIX_BINNING_H
