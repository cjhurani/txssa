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


#ifndef BLAS_FUNCTIONS_H
#define BLAS_FUNCTIONS_H

/* -------------------------------------------------------------------------- */

#include "fortran/fort_wrap_func.h"
#include "math/complex_types.h"

/* -------------------------------------------------------------------------- */

/* This header can be included in a C translation unit. */

/* Define the BLAS naming convention */
#define BLAS_double_type_id          d
#define BLAS_float_type_id           s
#define BLAS_complex_double_type_id  z
#define BLAS_complex_float_type_id   c

#define BLAS_double_TYPE_ID          D
#define BLAS_float_TYPE_ID           S
#define BLAS_complex_double_TYPE_ID  Z
#define BLAS_complex_float_TYPE_ID   C

typedef int BLAS_int;

/* Define the BLAS function declaration with types. */
#define BLAS_axpy_ARG(T)  BLAS_int* n, T* da, T* dx, BLAS_int* incx, T* dy, BLAS_int* incy
#define BLAS_copy_ARG(T)  BLAS_int* n, T* dx, BLAS_int* incx, T* dy, BLAS_int* incy
#define BLAS_gemm_ARG(T)  char* transa, char* transb, BLAS_int* m, BLAS_int* n, BLAS_int* k, T* alpha, T* a, BLAS_int* lda, T* b, BLAS_int* ldb, T* beta, T* c, BLAS_int* ldc
#define BLAS_gemv_ARG(T)  char* trans, BLAS_int* m, BLAS_int* n, T* alpha, T* a, BLAS_int* lda, T* x, BLAS_int* incx, T* beta, T* y, BLAS_int* incy
#define BLAS_syrk_ARG(T)  char* uplo, char* trans, BLAS_int* n, BLAS_int* k, T* alpha, T* a, BLAS_int* lda, T* beta, T* c, BLAS_int* ldc
#define BLAS_herk_ARG(T)  BLAS_syrk_ARG(T)
#define BLAS_symm_ARG(T)  char* side, char* uplo, BLAS_int* m, BLAS_int* n, T* alpha, T* a, BLAS_int* lda, T* b, BLAS_int* ldb, T* beta, T* c, BLAS_int* ldc
#define BLAS_trsm_ARG(T)  char* side, char* uplo, char* transa, char* diag, BLAS_int* m, BLAS_int* n, T* alpha, T* a, BLAS_int* lda, T* b, BLAS_int* ldb

/* -------------------------------------------------------------------------- */

/* Define the BLAS function declaration without types.  These will be used    */
/* for both the C and C++ function definitions.                               */
#define BLAS_copy_ARG_VAL  n, dx, incx, dy, incy
#define BLAS_axpy_ARG_VAL  n, da, dx, incx, dy, incy
#define BLAS_gemm_ARG_VAL  transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc
#define BLAS_gemv_ARG_VAL  trans, m, n, alpha, a, lda, x, incx, beta, y, incy
#define BLAS_syrk_ARG_VAL  uplo, trans, n, k, alpha, a, lda, beta, c, ldc
#define BLAS_herk_ARG_VAL  BLAS_syrk_ARG_VAL
#define BLAS_symm_ARG_VAL  side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc
#define BLAS_trsm_ARG_VAL  side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C"
{
#endif

/* Declare the C BLAS wrapper functions for all needed types. */
FORT_WRAP_FUNC_DECL_all(BLAS, axpy);
FORT_WRAP_FUNC_DECL_all(BLAS, copy);
FORT_WRAP_FUNC_DECL_all(BLAS, gemm);
FORT_WRAP_FUNC_DECL_all(BLAS, gemv);
FORT_WRAP_FUNC_DECL_all(BLAS, syrk);
FORT_WRAP_FUNC_DECL_complex(BLAS, herk);
FORT_WRAP_FUNC_DECL_all(BLAS, symm);
FORT_WRAP_FUNC_DECL_all(BLAS, trsm);

#ifdef __cplusplus
} // extern "C"
#endif

/* -------------------------------------------------------------------------- */

#endif /* BLAS_FUNCTIONS_H */
