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


#ifndef SPARSE_SPECTRAL_MISFIT_RHS_H
#define SPARSE_SPECTRAL_MISFIT_RHS_H

// -----------------------------------------------------------------------------

#include "sparse_vectors/sparse_vectors.h"
#include "math/precision_traits.h"
#include "math/complex_types.h"
#include "cpp/std_extensions.h"
#include "internal_api_error/internal_api_error.h"
#include <complex>
#include <cstddef>
#include <algorithm>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute the RHS of the derivative of quadratic form in misfit
// function J.
// -----------------------------------------------------------------------------

// See comments in sparse_spectral_misfit_lhs.h also.

// J = 0.5 * ( norm((X-A)*B1, 'fro')^2 + norm(B2*(X-A), 'fro')^2 )
// If X, A = m x n, B1 = n x q1, B2 = q2 x m, B2' * B2 = m x m, B1 * B1' = n x n.

// d1J(for X = 0) = (A * B1B1T + B2TB2 * A)      --> term-order is same as in J
// The above is true for any B1 and B2.  If B1 = B2 = pinv(A), then
// d1J(for X = 0) = 2 * pinv(A)', which is a special case.

// The code asks for the RHS matrix and fills the RHS vector for optimization.
// The caller needs to take care of the special cases and multiplying factors.
// The code works for a column matrix and column pattern.
// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type,
    typename extractor_type
>
bool sparse_spectral_misfit_rhs_internal(
    index_type num_rows,
    index_type num_cols,
    const value_type* RHS_col_values,  // num_rows x num_cols
    index_type RHS_col_leading_dim,
    offset_type num_dofs,
    const sparse_vectors_ids<index_type, offset_type>* col_split_pattern, // num_dofs, each num_rows vecs and num_cols size
    typename precision_traits<value_type>::scalar* b_values, // num_dofs
    const extractor_type& extractor)
{
    bool success =
        RHS_col_values &&
        num_rows <= RHS_col_leading_dim &&
        b_values &&
        col_split_pattern;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_rhs_internal: Unacceptable input argument(s).");

        return false;
    }

    typedef typename precision_traits<value_type>::scalar scalar_type;

    for(offset_type j_dof = 0; j_dof < num_dofs; ++j_dof)
    {
        const sparse_vectors_ids<index_type, offset_type>& j_dof_pat = col_split_pattern[j_dof];

        assert(j_dof_pat.num_vecs() == num_cols);
        assert(j_dof_pat.max_size() == num_rows);

        scalar_type tmp = 0;

        const value_type* A_off = RHS_col_values;

        for(index_type vec = 0; vec < num_cols; ++vec)
        {
            const index_type j_sz = j_dof_pat.num_vec_entries(vec);
            const index_type* j_ids = j_dof_pat.vec_ids_begin(vec);

            assert(j_sz <= num_cols);
            assert(j_ids);

            for(index_type j_id = 0; j_id < j_sz; ++j_id)
            {
                tmp += extractor(A_off[j_ids[j_id]]);
            }

            A_off += RHS_col_leading_dim;
        }

        b_values[j_dof] += tmp;
    }

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_spectral_misfit_rhs(
    index_type num_rows,
    index_type num_cols,
    const value_type* RHS_col_values,  // num_rows x num_cols
    index_type RHS_col_leading_dim,
    offset_type num_dofs,
    const sparse_vectors_ids<index_type, offset_type>* col_split_pattern, // num_dofs, each num_rows vecs and num_cols size
    value_type* b_values) // num_dofs
{
    bool success = b_values != 0;

    if(success)
    {
        std::fill(b_values, b_values + num_dofs, value_type(0));

        success =
            sparse_spectral_misfit_rhs_internal(
                num_rows, num_cols,
                RHS_col_values, RHS_col_leading_dim,
                num_dofs, col_split_pattern,
                b_values,
                std_extensions_identity<value_type>());
    }

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_rhs: Error.");
    }

    return success;
}

// -----------------------------------------------------------------------------

template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
bool sparse_spectral_misfit_rhs(
    index_type num_rows,
    index_type num_cols,
    const std::complex<scalar_type>* RHS_col_values,  // num_rows x num_cols
    index_type RHS_col_leading_dim,
    offset_type real_num_dofs,
    offset_type imag_num_dofs,
    const sparse_vectors_ids<index_type, offset_type>* real_col_split_pattern, // real_num_dofs, each num_rows vecs and num_cols size
    const sparse_vectors_ids<index_type, offset_type>* imag_col_split_pattern, // imag_num_dofs, each num_rows vecs and num_cols size
    scalar_type* b_values) // real_num_dofs + imag_num_dofs
{
    std::fill(b_values, b_values + std::size_t(real_num_dofs) + std::size_t(imag_num_dofs), scalar_type(0));

    return
        sparse_spectral_misfit_rhs_internal(
            num_rows, num_cols,
            RHS_col_values, RHS_col_leading_dim,
            real_num_dofs, real_col_split_pattern,
            b_values,
            real_extractor<scalar_type>())
        &&
        sparse_spectral_misfit_rhs_internal(
            num_rows, num_cols,
            RHS_col_values, RHS_col_leading_dim,
            imag_num_dofs, imag_col_split_pattern,
            b_values + real_num_dofs, // Imaginary part DOFs are after real.
            imag_extractor<scalar_type>());
}

// -----------------------------------------------------------------------------

#endif // SPARSE_SPECTRAL_MISFIT_RHS_H
