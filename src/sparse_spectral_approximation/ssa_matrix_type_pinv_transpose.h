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


#ifndef SSA_MATRIX_TYPE_PINV_TRANSPOSE_H
#define SSA_MATRIX_TYPE_PINV_TRANSPOSE_H

// -----------------------------------------------------------------------------

#include "txssa.h"
#include "dense_matrix_pinv/dense_matrix_qr_pinv.h"
#include "dense_vectors/dense_vectors.h"
#include "internal_api_error/internal_api_error.h"

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
bool ssa_matrix_type_pinv_transpose(
    index_type  num_rows,
    index_type  num_cols,
    value_type* A_col_values,
    index_type  A_col_leading_dim,
    ssa_matrix_type matrix_type,
    dense_vectors<index_type, value_type>* lnull,
    dense_vectors<index_type, value_type>* rnull)
{
    bool success =
        ssa_matrix_type_undefined < matrix_type &&
        matrix_type < ssa_matrix_type_num_types &&
        A_col_values &&
        num_rows <= A_col_leading_dim;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "ssa_matrix_type_pinv_transpose: Unacceptable input argument(s).");

        return false;
    }

    //TODO: Use type dependent pinv algorithm.

    // if(matrix_type == ssa_matrix_type_general)
    {
        success = dense_matrix_qr_pinv_transpose(
            num_rows, num_cols,
            A_col_values, A_col_leading_dim,
            lnull, rnull);
    }

    if(!success)
        internal_api_error_set_last(
            "ssa_matrix_type_pinv_transpose: Error.");

    return success;
}

// -----------------------------------------------------------------------------

#endif // SSA_MATRIX_TYPE_PINV_TRANSPOSE_H
