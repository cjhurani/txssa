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


#ifndef SPARSE_SPECTRAL_MISFIT_LHS_MATRICES_H
#define SPARSE_SPECTRAL_MISFIT_LHS_MATRICES_H

// -----------------------------------------------------------------------------

#include "dense_algorithms/dense_matrix_utils.h"
#include "blas_wrap/dense_matrix_mult.h"
#include "math/complex_types.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute one or both matrix products for the LHS of the
// minimization problem (A*A' and A'*A).  Even though the output is
// Hermitian, we still fill both halves.
// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool sparse_spectral_misfit_lhs_matrices(
// Input:
    index_type num_rows,
    index_type num_cols,
    const value_type* A_col_values, // num_rows x num_cols
    index_type A_col_leading_dim,

// Output:
    value_type* AAT_col_values,  // num_rows x num_rows
    index_type  AAT_col_leading_dim,
    value_type* ATA_col_values,  // num_cols x num_cols
    index_type  ATA_col_leading_dim)
{
    bool success =
        A_col_values &&
        num_rows <= A_col_leading_dim &&
        (AAT_col_values  || ATA_col_values) &&
        (!AAT_col_values || num_rows <= AAT_col_leading_dim) &&
        (!ATA_col_values || num_cols <= ATA_col_leading_dim);

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_lhs_matrices: Unacceptable input argument(s).");

        return false;
    }

    value_type (*complex_conjugate)(const value_type&) = std::conj;

    if(AAT_col_values)
    {
        success =
            dense_matrix_mult_herk(
                'U', 'N',
                num_rows, num_cols,
                value_type(1),
                A_col_values, A_col_leading_dim,
                value_type(0),
                AAT_col_values, AAT_col_leading_dim)
            &&
            dense_matrix_utils_copy_upper_to_lower_func( // Only 'U' above
                num_rows,
                AAT_col_values, AAT_col_leading_dim,
                complex_conjugate);
    }

    if(ATA_col_values)
    {
        success =
            dense_matrix_mult_herk(
                'U', 'C',
                num_cols, num_rows,
                value_type(1),
                A_col_values, A_col_leading_dim,
                value_type(0),
                ATA_col_values, ATA_col_leading_dim)
            &&
            dense_matrix_utils_copy_upper_to_lower_func( // Only 'U' above
                num_cols,
                ATA_col_values, ATA_col_leading_dim,
                complex_conjugate);
    }

    if(!success)
    {
        assert(false);
        internal_api_error_set_last(
            "sparse_spectral_misfit_lhs_matrices: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

#endif // SPARSE_SPECTRAL_MISFIT_LHS_MATRICES_H
