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


#ifndef LAPACK_FUNCTIONS_H
#define LAPACK_FUNCTIONS_H

/* -------------------------------------------------------------------------- */

#include "fortran/fort_wrap_func.h"
#include "math/complex_types.h"

/* -------------------------------------------------------------------------- */

/* This header can be included in a C translation unit. */

/* Define the LAPACK naming convention */
#define LAPACK_double_type_id          d
#define LAPACK_float_type_id           s
#define LAPACK_complex_double_type_id  z
#define LAPACK_complex_float_type_id   c

#define LAPACK_double_TYPE_ID          D
#define LAPACK_float_TYPE_ID           S
#define LAPACK_complex_double_TYPE_ID  Z
#define LAPACK_complex_float_TYPE_ID   C

/* Use for numbers only, not for info argument. */
typedef int LAPACK_int;

/* Define the LAPACK function declaration with types. */
#define LAPACK_trtri_ARG(T)  char* uplo, char* diag, LAPACK_int* n, T* a, LAPACK_int* lda, int* info
#define LAPACK_geqrf_ARG(T)  LAPACK_int* m, LAPACK_int* n, T* a, LAPACK_int* lda, T* tau, T* work, LAPACK_int* lwork, int* info
#define LAPACK_potrf_ARG(T)  char* uplo, LAPACK_int* n, T* a, LAPACK_int* lda, int* info
#define LAPACK_posv_ARG(T)   char* uplo, LAPACK_int* n, LAPACK_int* nrhs, T* a, LAPACK_int* lda, T* b, LAPACK_int* ldb, int* info
#define LAPACK_orgqr_ARG(T)  LAPACK_int* m, LAPACK_int* n, LAPACK_int* k, T* a, LAPACK_int* lda, T* tau, T* work, LAPACK_int* lwork, int* info
#define LAPACK_ungqr_ARG(T)  LAPACK_orgqr_ARG(T)
#define LAPACK_ormqr_ARG(T)  char* side, char* trans, LAPACK_int* m, LAPACK_int* n, LAPACK_int* k, T* a, LAPACK_int* lda, T* tau, T* c, LAPACK_int* ldc, T* work, LAPACK_int* lwork, int* info
#define LAPACK_unmqr_ARG(T)  LAPACK_ormqr_ARG(T)

#define LAPACK_geqp3_ARG(T)           LAPACK_int* m, LAPACK_int* n, T* a, LAPACK_int* lda, LAPACK_int* jpvt, T* tau, T* work, LAPACK_int* lwork, int* info
#define LAPACK_geqp3_ARG_2(T, T_REAL) LAPACK_int* m, LAPACK_int* n, T* a, LAPACK_int* lda, LAPACK_int* jpvt, T* tau, T* work, LAPACK_int* lwork, T_REAL* rwork, int* info

#define LAPACK_gesvd_ARG(T)           char* jobu, char* jobvt, LAPACK_int* m, LAPACK_int* n, T* a, LAPACK_int* lda, T* s, T* u, LAPACK_int* ldu, T* vt, LAPACK_int* ldvt, T* work, LAPACK_int* lwork, int* info
#define LAPACK_gesvd_ARG_2(T, T_REAL) char* jobu, char* jobvt, LAPACK_int* m, LAPACK_int* n, T* a, LAPACK_int* lda, T_REAL* s, T* u, LAPACK_int* ldu, T* vt, LAPACK_int* ldvt, T* work, LAPACK_int* lwork, T_REAL* rwork, int* info


/* -------------------------------------------------------------------------- */

/* Define the LAPACK function declaration without types.  These will be used  */
/* for both the C and C++ function definitions.                               */
#define LAPACK_trtri_ARG_VAL  uplo, diag, n, a, lda, info
#define LAPACK_geqrf_ARG_VAL  m, n, a, lda, tau, work, lwork, info
#define LAPACK_potrf_ARG_VAL  uplo, n, a, lda, info
#define LAPACK_posv_ARG_VAL   uplo, n, nrhs, a, lda, b, ldb, info
#define LAPACK_orgqr_ARG_VAL  m, n, k, a, lda, tau, work, lwork, info
#define LAPACK_ungqr_ARG_VAL  LAPACK_orgqr_ARG_VAL
#define LAPACK_ormqr_ARG_VAL  side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info
#define LAPACK_unmqr_ARG_VAL  LAPACK_ormqr_ARG_VAL

#define LAPACK_geqp3_ARG_VAL    m, n, a, lda, jpvt, tau, work, lwork, info
#define LAPACK_geqp3_ARG_2_VAL  m, n, a, lda, jpvt, tau, work, lwork, rwork, info

#define LAPACK_gesvd_ARG_VAL    jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info
#define LAPACK_gesvd_ARG_2_VAL  jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C"
{
#endif

/* Declare the C LAPACK wrapper functions for all needed types. */
FORT_WRAP_FUNC_DECL_all(          LAPACK, trtri);
FORT_WRAP_FUNC_DECL_all(          LAPACK, geqrf);
FORT_WRAP_FUNC_DECL_all(          LAPACK, potrf);
FORT_WRAP_FUNC_DECL_all(          LAPACK, posv);
FORT_WRAP_FUNC_DECL_real(         LAPACK, orgqr);
FORT_WRAP_FUNC_DECL_complex(      LAPACK, ungqr);
FORT_WRAP_FUNC_DECL_real(         LAPACK, ormqr);
FORT_WRAP_FUNC_DECL_complex(      LAPACK, unmqr);
FORT_WRAP_FUNC_DECL_real(         LAPACK, geqp3);
FORT_WRAP_FUNC_DECL_ARG_2_complex(LAPACK, geqp3);
FORT_WRAP_FUNC_DECL_real(         LAPACK, gesvd);
FORT_WRAP_FUNC_DECL_ARG_2_complex(LAPACK, gesvd);

#ifdef __cplusplus
} // extern "C"
#endif

/* -------------------------------------------------------------------------- */

#endif /* LAPACK_FUNCTIONS_H */
