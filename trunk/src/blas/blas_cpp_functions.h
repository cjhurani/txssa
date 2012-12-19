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


#ifndef BLAS_CPP_FUNCTIONS_H
#define BLAS_CPP_FUNCTIONS_H

// -----------------------------------------------------------------------------

#ifdef __cplusplus

// -----------------------------------------------------------------------------

#include "blas/blas_functions.h"
#include "fortran/fort_wrap_cpp_func.h"

// -----------------------------------------------------------------------------

// Define the C++ BLAS function definition so that overloading works.
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, axpy)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, copy)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, gemm)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, gemv)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, syrk)
FORT_WRAP_CPP_FUNC_DEF_complex(BLAS, herk)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, symm)
FORT_WRAP_CPP_FUNC_DEF_all(BLAS, trsm)

// -----------------------------------------------------------------------------

// BLAS uses "sy" for real-symmetric and complex-symmetric and "he" for
// complex-Hermitian.  Mathematically and programmatically it leads to problems
// because it is better to have one function name for real-symmetric and
// complex-Hermitian and a different function name for complex-symmetric.
// Thus, we create and overload "he" for the real case so that it works for
// real-symmetric.  To avoid surprises, we don't do this at the C wrapper level 
// and because function overloading wouldn't work anyway.  It is done at C++
// level.  The idea is to call BLAS_herk/BLAS_hemm/etc for real as well as
// complex matrices if the matrix is Hermitian (a generality of symmetric).

// inline for weak symbols

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(BLAS, herk)(BLAS_syrk_ARG(double))
{
    FORT_WRAP_FUNC(BLAS, double, syrk)(BLAS_syrk_ARG_VAL);
}

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(BLAS, herk)(BLAS_syrk_ARG(float))
{
    FORT_WRAP_FUNC(BLAS, float, syrk)(BLAS_syrk_ARG_VAL);
}

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(BLAS, hemm)(BLAS_symm_ARG(double))
{
    FORT_WRAP_FUNC(BLAS, double, symm)(BLAS_symm_ARG_VAL);
}

inline FORT_RET FORT_WRAP_FUNC_NOTYPE(BLAS, hemm)(BLAS_symm_ARG(float))
{
    FORT_WRAP_FUNC(BLAS, float, symm)(BLAS_symm_ARG_VAL);
}

// -----------------------------------------------------------------------------

#endif // __cplusplus

// -----------------------------------------------------------------------------

#endif // BLAS_CPP_FUNCTIONS_H
