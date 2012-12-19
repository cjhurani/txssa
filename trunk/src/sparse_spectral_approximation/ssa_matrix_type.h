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


#ifndef SSA_MATRIX_TYPE_H
#define SSA_MATRIX_TYPE_H

/* -------------------------------------------------------------------------- */

#include "txssa.h"

/* This header can be included in a C translation unit. */

#ifdef __cplusplus
extern "C" {
#endif

int ssa_matrix_type_is_hermitian                (enum ssa_matrix_type type);
int ssa_matrix_type_is_abs_sym                  (enum ssa_matrix_type type);
int ssa_matrix_type_is_left_right_nullity_equal (enum ssa_matrix_type type);
int ssa_matrix_type_is_normal                   (enum ssa_matrix_type type);
int ssa_matrix_type_is_AAT_computable_from_ATA  (enum ssa_matrix_type type);
int ssa_matrix_type_is_real_part_symmetric      (enum ssa_matrix_type type);
int ssa_matrix_type_is_imag_part_symmetric      (enum ssa_matrix_type type);

#ifdef __cplusplus
} // extern "C"
#endif

/* -------------------------------------------------------------------------- */

#endif /* SSA_MATRIX_TYPE_H */
