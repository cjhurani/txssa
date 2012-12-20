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


#ifndef LAPACK_CPP_FUNCTIONS_H
#define LAPACK_CPP_FUNCTIONS_H

// -----------------------------------------------------------------------------

#ifdef __cplusplus

// -----------------------------------------------------------------------------

#include "lapack/lapack_functions.h"
#include "fortran/fort_wrap_cpp_func.h"

// -----------------------------------------------------------------------------

// Define the C++ LAPACK function definition so that overloading works.
FORT_WRAP_CPP_FUNC_DEF_all(LAPACK, trtri)
FORT_WRAP_CPP_FUNC_DEF_all(LAPACK, geqrf)
FORT_WRAP_CPP_FUNC_DEF_all(LAPACK, potrf)
FORT_WRAP_CPP_FUNC_DEF_all(LAPACK, posv)

// LAPACK has [DS]ORGQR (orthogonal) but the corresponding complex subroutines
// are called [ZC]UNGQR (unitary).  For consistency in the C++ version so that
// a single LAPACK_ungqr works for all types, we overload LAPACK_ungqr.
// Same logic for ormqr/unmqr.

FORT_WRAP_CPP_FUNC_DEF_real(LAPACK, orgqr)
FORT_WRAP_CPP_FUNC_DEF_complex(LAPACK, ungqr)

// inline for weak symbols
inline FORT_RET FORT_WRAP_FUNC_NOTYPE(LAPACK, ungqr)(LAPACK_ungqr_ARG(double))
{
    FORT_WRAP_FUNC(LAPACK, double, orgqr)(LAPACK_orgqr_ARG_VAL);
}

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(LAPACK, ungqr)(LAPACK_ungqr_ARG(float))
{
    FORT_WRAP_FUNC(LAPACK, float,  orgqr)(LAPACK_orgqr_ARG_VAL);
}

FORT_WRAP_CPP_FUNC_DEF_real(LAPACK, ormqr)
FORT_WRAP_CPP_FUNC_DEF_complex(LAPACK, unmqr)

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(LAPACK, unmqr)(LAPACK_unmqr_ARG(double))
{
    FORT_WRAP_FUNC(LAPACK, double, ormqr)(LAPACK_ormqr_ARG_VAL);
}

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(LAPACK, unmqr)(LAPACK_unmqr_ARG(float))
{
    FORT_WRAP_FUNC(LAPACK, float,  ormqr)(LAPACK_ormqr_ARG_VAL);
}

// LAPACK takes extra arguments of a different type for [ZC]GEQP3 compared to
// the [DS] versions.  For consistency in the C++ version so that a single
// LAPACK_geqp3 works for all types, we overload LAPACK_geqp3 to take extra
// arguments for the [DS] version and ignore them.  One can pass null pointers
// for the extra arguments in [DS] cases.

FORT_WRAP_CPP_FUNC_DEF_real(LAPACK, geqp3)
FORT_WRAP_CPP_FUNC_DEF_ARG_2_complex(LAPACK, geqp3)

inline FORT_RET
FORT_WRAP_FUNC_NOTYPE(LAPACK, geqp3)(LAPACK_geqp3_ARG_2(double, double))
{
    (void) rwork;
    FORT_WRAP_FUNC(LAPACK, double, geqp3)(LAPACK_geqp3_ARG_VAL);
}

inline FORT_RET
FORT_WRAP_FUNC_NOTYPE(LAPACK, geqp3)(LAPACK_geqp3_ARG_2(float,  float))
{
    (void) rwork;
    FORT_WRAP_FUNC(LAPACK, float,  geqp3)(LAPACK_geqp3_ARG_VAL);
}

// -----------------------------------------------------------------------------

FORT_WRAP_CPP_FUNC_DEF_real(LAPACK, gesvd)
FORT_WRAP_CPP_FUNC_DEF_ARG_2_complex(LAPACK, gesvd)

inline FORT_RET
FORT_WRAP_FUNC_NOTYPE(LAPACK, gesvd)(LAPACK_gesvd_ARG_2(double, double))
{
    (void) rwork;
    FORT_WRAP_FUNC(LAPACK, double, gesvd)(LAPACK_gesvd_ARG_VAL);
}

inline FORT_RET
FORT_WRAP_FUNC_NOTYPE(LAPACK, gesvd)(LAPACK_gesvd_ARG_2(float,  float))
{
    (void) rwork;
    FORT_WRAP_FUNC(LAPACK, float,  gesvd)(LAPACK_gesvd_ARG_VAL);
}

// -----------------------------------------------------------------------------

#endif // __cplusplus

// -----------------------------------------------------------------------------

#endif // LAPACK_CPP_FUNCTIONS_H
