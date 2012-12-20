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


#ifndef SPARSE_SPECTRAL_MISFIT_LHS_H
#define SPARSE_SPECTRAL_MISFIT_LHS_H

// -----------------------------------------------------------------------------

#include "sparse_vectors/sparse_vectors.h"
#include "dense_algorithms/dense_matrix_utils.h"
#include "math/precision_traits.h"
#include "math/complex_types.h"
#include "cpp/std_extensions.h"
#include "internal_api_error/internal_api_error.h"
#include <complex>
#include <cstddef>
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Compute the Hessian of the quadratic form in misfit function J.
// -----------------------------------------------------------------------------

// J = 0.5 * ( norm((X-A)*B1, 'fro')^2 + norm(B2*(X-A), 'fro')^2 )
// If X, A = m x n, B1 = n x q1, B2 = q2 x m, B2' * B2 = m x m, B1 * B1' = n x n.

// d2J = kron(pinv_ATA, eye(m)) + kron(eye(n), pinv_AAT)  (if B1 = B2 = pinv(A))
// pinv_ATA = pinv(A'A) = pinv(A)pinv(A') = B1B1T
// pinv_AAT = pinv(AA') = pinv(A')pinv(A) = B2TB2
// d2J = kron(B1B1T, eye(m)) + kron(eye(n), B2TB2)  --> term-order is same as J

// -----------------------------------------------------------------------------

// This function, when called appropriately, can work for real and
// complex and also different Frobenius norms.

template
<
    typename index_type,
    typename offset_type,
    typename value_type,
    typename extractor_type
>
bool sparse_spectral_misfit_lhs_internal(
    index_type num_vecs,
    index_type max_size,
    const value_type* quad_col_values,  // max_size x max_size, must be full.
    index_type quad_col_leading_dim,
    offset_type num_dofs_1,
    offset_type num_dofs_2,
    const sparse_vectors_ids<index_type, offset_type>* split_pat_1, // num_dofs_1, each num_vecs vecs and max_size size
    const sparse_vectors_ids<index_type, offset_type>* split_pat_2, // num_dofs_2, each num_vecs vecs and max_size size
    typename precision_traits<value_type>::scalar* LS_A_col_values, // num_dofs_1 x num_dofs_2
    offset_type LS_A_col_leading_dim,
    const extractor_type& extractor,
    bool upper_half_only = true)
{
    bool success =
        quad_col_values &&
        LS_A_col_values &&
        max_size <= quad_col_leading_dim &&
        num_dofs_1 <= LS_A_col_leading_dim &&
        split_pat_1 &&
        split_pat_2 &&
        ((num_dofs_1 == num_dofs_2) || !upper_half_only); // Need to ask for full if LS_A is not square (just an extra check)

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_lhs_internal: Unacceptable input argument(s).");

        return false;
    }

    typedef typename precision_traits<value_type>::scalar scalar_type;

    // Go through columns of LS_A_col_values (which is to be updated).
    // Each entry in (upper half or full) LS_A_col_values will be updated once only.

    for(offset_type j_dof = 0; j_dof < num_dofs_2; ++j_dof)
    {
        const sparse_vectors_ids<index_type, offset_type>& j_dof_pat = split_pat_2[j_dof];

        assert(j_dof_pat.num_vecs() == num_vecs);
        assert(j_dof_pat.max_size() == max_size);

        const offset_type i_dof_end = upper_half_only ? j_dof + 1 : num_dofs_1;

        for(offset_type i_dof = 0; i_dof < i_dof_end; ++i_dof) // upper triangle or full LS_A
        {
            const sparse_vectors_ids<index_type, offset_type>& i_dof_pat = split_pat_1[i_dof];

            assert(i_dof_pat.num_vecs() == num_vecs);
            assert(i_dof_pat.max_size() == max_size);

            scalar_type tmp = 0;

            for(index_type vec = 0; vec < num_vecs; ++vec)
            {
                const index_type j_sz = j_dof_pat.num_vec_entries(vec);

                if(j_sz)
                {
                    const index_type i_sz = i_dof_pat.num_vec_entries(vec);

                    if(i_sz)
                    {
                        const index_type* j_ids = j_dof_pat.vec_ids_begin(vec);
                        const index_type* i_ids = i_dof_pat.vec_ids_begin(vec);

                        for(index_type j_id = 0; j_id < j_sz; ++j_id)
                        {
                            const value_type* quad_off = quad_col_values +
                                std::size_t(j_ids[j_id]) * std::size_t(quad_col_leading_dim);

                            for(index_type i_id = 0; i_id < i_sz; ++i_id)
                            {
                                tmp += extractor(quad_off[i_ids[i_id]]);
                            }
                        }
                    }
                }
            }

            LS_A_col_values[i_dof] += tmp;
        }

        LS_A_col_values += LS_A_col_leading_dim;
    }

    return success;
}

// -----------------------------------------------------------------------------

// B1B1T_col_values or B2TB2_col_values can be 0 but not both.
// For the 0 case, the corresponding col_leading_dim is not used.

template
<
    typename index_type,
    typename offset_type,
    typename value_type
>
bool sparse_spectral_misfit_lhs(
// Input:
    index_type num_rows,
    index_type num_cols,
    const value_type* B2TB2_col_values,  // num_rows x num_rows
    index_type  B2TB2_col_leading_dim,
    const value_type* B1B1T_col_values,  // num_cols x num_cols
    index_type  B1B1T_col_leading_dim,
    offset_type num_dofs,
    const sparse_vectors_ids<index_type, offset_type>* row_split_pattern, // num_dofs
    const sparse_vectors_ids<index_type, offset_type>* col_split_pattern, // num_dofs

// Output:
    value_type* LS_A_col_values, // num_dofs x num_dofs
    offset_type LS_A_col_leading_dim)
{
    // Only upper triangle of LS_A_col_values will be filled.

    bool success =
        B2TB2_col_values || B1B1T_col_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_lhs: Unacceptable input argument(s).");

        return false;
    }

    return
        dense_matrix_utils_fill_upper(
            num_dofs, num_dofs,
            LS_A_col_values, LS_A_col_leading_dim,
            value_type(0))
        &&
        (B1B1T_col_values ? sparse_spectral_misfit_lhs_internal(
            num_rows, num_cols,
            B1B1T_col_values, B1B1T_col_leading_dim,
            num_dofs, num_dofs, row_split_pattern, row_split_pattern,
            LS_A_col_values, LS_A_col_leading_dim,
            std_extensions_identity<value_type>()) : true)
        &&
        (B2TB2_col_values ? sparse_spectral_misfit_lhs_internal(
            num_cols, num_rows,
            B2TB2_col_values, B2TB2_col_leading_dim,
            num_dofs, num_dofs, col_split_pattern, col_split_pattern,
            LS_A_col_values, LS_A_col_leading_dim,
            std_extensions_identity<value_type>()) : true);
}

// -----------------------------------------------------------------------------

// B1B1T_col_values or B2TB2_col_values can be 0 but not both.
// For the 0 case, the corresponding col_leading_dim is not used.

template
<
    typename index_type,
    typename offset_type,
    typename scalar_type
>
bool sparse_spectral_misfit_lhs(
// Input:
    index_type num_rows,
    index_type num_cols,
    const std::complex<scalar_type>* B2TB2_col_values,  // num_rows x num_rows
    index_type  B2TB2_col_leading_dim,
    const std::complex<scalar_type>* B1B1T_col_values,  // num_cols x num_cols
    index_type  B1B1T_col_leading_dim,
    offset_type real_num_dofs,
    offset_type imag_num_dofs,
    const sparse_vectors_ids<index_type, offset_type>* real_row_split_pattern, // real_num_dofs
    const sparse_vectors_ids<index_type, offset_type>* imag_row_split_pattern, // imag_num_dofs
    const sparse_vectors_ids<index_type, offset_type>* real_col_split_pattern, // real_num_dofs
    const sparse_vectors_ids<index_type, offset_type>* imag_col_split_pattern, // imag_num_dofs

// Output:
    scalar_type* LS_A_col_values, // (real_num_dofs + imag_num_dofs)^2
    offset_type LS_A_col_leading_dim)
{
    // Only upper triangle of LS_A_col_values will be filled.

    bool success =
        B2TB2_col_values || B1B1T_col_values;

    if(!success)
    {
        assert(false);

        internal_api_error_set_last(
            "sparse_spectral_misfit_lhs: Unacceptable input argument(s).");

        return false;
    }

    const offset_type num_dofs = real_num_dofs + imag_num_dofs;

    // For each B (1 or 2), we need two halves (real <-> real and
    // imag <-> imag) and one full matrix for upper half (real <-> imag).
    // Hence total 6 calls to the template routine.

    const std::size_t offset_01 = std::size_t(LS_A_col_leading_dim) * std::size_t(real_num_dofs);
    const std::size_t offset_11 = offset_01 + std::size_t(real_num_dofs);

    return
        dense_matrix_utils_fill_upper(
            num_dofs, num_dofs,
            LS_A_col_values, LS_A_col_leading_dim,
            scalar_type(0))
        &&
        (B1B1T_col_values ?
            sparse_spectral_misfit_lhs_internal(
                num_rows, num_cols,
                B1B1T_col_values, B1B1T_col_leading_dim,
                real_num_dofs, real_num_dofs, real_row_split_pattern, real_row_split_pattern,
                LS_A_col_values, LS_A_col_leading_dim,
                real_extractor<scalar_type>())
            &&
            sparse_spectral_misfit_lhs_internal(
                num_rows, num_cols,
                B1B1T_col_values, B1B1T_col_leading_dim,
                imag_num_dofs, imag_num_dofs, imag_row_split_pattern, imag_row_split_pattern,
                LS_A_col_values + offset_11, LS_A_col_leading_dim,
                real_extractor<scalar_type>())
            &&
            sparse_spectral_misfit_lhs_internal(
                num_rows, num_cols,
                B1B1T_col_values, B1B1T_col_leading_dim,
                real_num_dofs, imag_num_dofs, real_row_split_pattern, imag_row_split_pattern,
                LS_A_col_values + offset_01, LS_A_col_leading_dim,
                imag_extractor<scalar_type>(), false) // full
            : true )
        &&
        (B2TB2_col_values ?
            sparse_spectral_misfit_lhs_internal(
                num_cols, num_rows,
                B2TB2_col_values, B2TB2_col_leading_dim,
                real_num_dofs, real_num_dofs, real_col_split_pattern, real_col_split_pattern,
                LS_A_col_values, LS_A_col_leading_dim,
                real_extractor<scalar_type>())
            &&
            sparse_spectral_misfit_lhs_internal(
                num_cols, num_rows,
                B2TB2_col_values, B2TB2_col_leading_dim,
                imag_num_dofs, imag_num_dofs, imag_col_split_pattern, imag_col_split_pattern,
                LS_A_col_values + offset_11, LS_A_col_leading_dim,
                real_extractor<scalar_type>())
            &&
            sparse_spectral_misfit_lhs_internal(
                num_cols, num_rows,
                B2TB2_col_values, B2TB2_col_leading_dim,
                real_num_dofs, imag_num_dofs, real_col_split_pattern, imag_col_split_pattern,
                LS_A_col_values + offset_01, LS_A_col_leading_dim,
                conj_imag_extractor<scalar_type>(), false) // full
            : true );
}

// -----------------------------------------------------------------------------

#endif // SPARSE_SPECTRAL_MISFIT_LHS_H
