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


// -----------------------------------------------------------------------------

#include "txssa.h"
#include "sparse_spectral_approximation/ssa_matrix_type.h"
#include "internal_api_error/internal_api_error.h"
#include <cassert>

// -----------------------------------------------------------------------------

#ifdef _MSC_VER
#pragma warning( disable : 4514 ) // unreferenced inline function has been removed
#pragma warning( disable : 4710 ) // function not inlined
#endif

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Objective: Provide information on whether a matrix type has a given property.
// -----------------------------------------------------------------------------

// anonymous
namespace
{
    enum matrix_property
    {
        // Do not change the order unless you also change
        // matrix_type_to_properties array below.

        matrix_property_undefined = -1,
        matrix_property_hermitian,
        matrix_property_abs_sym,
        matrix_property_left_right_nullity_equal,
        matrix_property_normal,
        matrix_property_AAT_computable_from_ATA,
        matrix_property_real_part_symmetric,
        matrix_property_imag_part_symmetric,
        matrix_property_num_properties
    };

    const bool matrix_type_to_properties[ssa_matrix_type_num_types][matrix_property_num_properties] =
    {
        {0,0,0,0,0,0,0}, // general
        {1,1,1,1,1,1,0}, // hermitian_pos_def
        {1,1,1,1,1,1,0}, // hermitian_pos_semi_def
        {1,1,1,1,1,1,0}, // hermitian
        {0,1,1,1,1,0,1}, // skew_hermitian
        {0,1,1,0,1,1,1}  // complex_symmetric
    };

    // -----------------------------------------------------------------------------

    int ssa_matrix_type_has_property(
        ssa_matrix_type type,
        matrix_property prop)
    {
        bool success =
            ssa_matrix_type_undefined < type &&
            type < ssa_matrix_type_num_types &&
            matrix_property_undefined < prop &&
            prop < matrix_property_num_properties;

        int ret = 0;

        if(success)
        {
            ret = matrix_type_to_properties[type][prop];
        }
        else
        {
            assert(false);
            internal_api_error_set_last(
                "ssa_matrix_type_has_property: Unacceptable input argument(s).");
        }

        return ret;
    }

} // namespace

// -----------------------------------------------------------------------------

int ssa_matrix_type_is_hermitian(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_hermitian);
}

int ssa_matrix_type_is_abs_sym(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_abs_sym);
}

int ssa_matrix_type_is_left_right_nullity_equal(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_left_right_nullity_equal);
}

int ssa_matrix_type_is_normal(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_normal);
}

int ssa_matrix_type_is_AAT_computable_from_ATA(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_AAT_computable_from_ATA);
}

int ssa_matrix_type_is_real_part_symmetric(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_real_part_symmetric);
}

int ssa_matrix_type_is_imag_part_symmetric(ssa_matrix_type type)
{
    return ssa_matrix_type_has_property(type, matrix_property_imag_part_symmetric);
}

// -----------------------------------------------------------------------------
