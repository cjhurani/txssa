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


#ifndef SPLIT_PATTERN_TO_BINS_H
#define SPLIT_PATTERN_TO_BINS_H

// -----------------------------------------------------------------------------

#include "sparse_vectors/sparse_vectors.h"
#include "internal_api_error/internal_api_error.h"
#include <vector>
#include <stdexcept>
#include <algorithm> // std::fill
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Split a matrix with integral values into multiple patterns
// based on the values.
// -----------------------------------------------------------------------------

template<typename index_type, typename offset_type>
bool split_pattern_to_bins(
    index_type   n_vecs,
    index_type   max_vec_size,
    const offset_type* offsets,
    const index_type*  ids,
    const offset_type* values,
    offset_type  num_bins,
    sparse_vectors_ids<index_type, offset_type>* split) // num_bins size
{
    bool success =
        offsets &&
        ((!n_vecs || !max_vec_size) || values) &&
        split;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "split_pattern_to_bins: Unacceptable input argument(s).");

        return false;
    }

    std::vector< std::vector<index_type> > counts_per_bin_per_vec;

    try
    {
        counts_per_bin_per_vec.resize(num_bins, std::vector<index_type>(n_vecs));
    }
    catch(const std::exception& exc)
    {
        assert(false);
        internal_api_error_set_last(
            (std::string("split_pattern_to_bins: Exception. ") + exc.what()));
        return false;
    }

    // Go through the pattern to compute counts_per_bin_per_vec

    for(index_type i = 0; i < n_vecs; ++i)
    {
        for(offset_type j = offsets[i]; j < offsets[i+1]; ++j)
        {
            const offset_type bin_id = values[j];
            assert(bin_id < num_bins);

            ++counts_per_bin_per_vec[bin_id][i];
        }
    }

    success = true;

    for(offset_type k = 0; (k < num_bins) && success; ++k)
    {
        success = split[k].allocate(n_vecs, max_vec_size,
            &counts_per_bin_per_vec[k].front());
        assert(success);
    }

    if(success)
    {
        // Reuse the counts_per_bin_per_vec space.
        // First fill with 0 and then use it to track the "next" location
        // to be filled in the loop below.

        for(offset_type bin_id = 0; bin_id < num_bins; ++bin_id)
        {
            std::fill(
                counts_per_bin_per_vec[bin_id].begin(),
                counts_per_bin_per_vec[bin_id].end(),
                0);
        }

        // Go through the pattern again, and fill the
        // offsets in each vec and each bin.

        for(index_type i = 0; i < n_vecs; ++i)
        {
            for(offset_type j = offsets[i]; j < offsets[i+1]; ++j)
            {
                const offset_type bin_id = values[j];

                index_type& next_for_this_bin = counts_per_bin_per_vec[bin_id][i];

                assert(next_for_this_bin < split[bin_id].num_vec_entries(i));

                assert(ids[j] < max_vec_size);

                split[bin_id].vec_ids_begin(i)[next_for_this_bin++] = ids[j];
            }
        }
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "split_pattern_to_bins: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPLIT_PATTERN_TO_BINS_H
