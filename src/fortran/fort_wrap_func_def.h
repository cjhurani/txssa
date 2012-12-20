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


/*
 This file should be included by the implementing translation unit of the
 wrapper code.  For this reason this file allows being included only once.

 This file requires FORTRAN_SYMBOL_SCHEME to be defined appropriately
 as an integer to compile successfully.
 */

/* -------------------------------------------------------------------------- */

#include "fortran/fort_wrap_func.h"

/* -------------------------------------------------------------------------- */

#ifdef FORT_FUNC_NAME

#error "FORT_FUNC_NAME already defined.  Note: this file must NOT be included more than once in a translation unit."

#endif

/* -------------------------------------------------------------------------- */

#if defined(FORTRAN_SYMBOL_SCHEME) && (1 <= FORTRAN_SYMBOL_SCHEME) && (FORTRAN_SYMBOL_SCHEME <= 4)


#if   FORTRAN_SYMBOL_SCHEME == 1

#define FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID) type_id ## name ## _

#elif FORTRAN_SYMBOL_SCHEME == 2

#define FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID) type_id ## name

#elif FORTRAN_SYMBOL_SCHEME == 3

#define FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID) TYPE_ID ## NAME ## _

#elif FORTRAN_SYMBOL_SCHEME == 4

#define FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID) TYPE_ID ## NAME

// More schemes can be put here.

#else

#error "Internal error - should never reach here."

#endif


#else

#error "Either FORTRAN_SYMBOL_SCHEME is not defined, or it is out of valid range. dcopy_=1, dcopy=2, DCOPY_=3, DCOPY=4."

#endif /* FORTRAN_SYMBOL_SCHEME */

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus

#define FORT_WRAP_EXTERN_C_BEGIN extern "C" {
#define FORT_WRAP_EXTERN_C_END   }

#else

#define FORT_WRAP_EXTERN_C_BEGIN
#define FORT_WRAP_EXTERN_C_END

#endif

/* -------------------------------------------------------------------------- */

#ifdef FORT_FUNC_NAME

#define FORT_WRAP_FUNC_decl_and_def(lib, type, name, NAME, type_id, TYPE_ID)                          \
    FORT_WRAP_EXTERN_C_BEGIN                                                                          \
        FORT_RET FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID)(FORT_WRAP_FUNC_ARG(lib, type, name));   \
        FORT_RET FORT_WRAP_FUNC(lib, type, name)(FORT_WRAP_FUNC_ARG(lib, type, name))                 \
        {                                                                                             \
            FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID)(FORT_WRAP_FUNC_ARG_VAL(lib, name));          \
        }                                                                                             \
    FORT_WRAP_EXTERN_C_END

#define FORT_WRAP_FUNC_decl_and_def_2(lib, type, type_2, name, NAME, type_id, TYPE_ID)                         \
    FORT_WRAP_EXTERN_C_BEGIN                                                                                   \
        FORT_RET FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID)(FORT_WRAP_FUNC_ARG_2(lib, type, type_2, name));  \
        FORT_RET FORT_WRAP_FUNC(lib, type, name)(FORT_WRAP_FUNC_ARG_2(lib, type, type_2, name))                \
        {                                                                                                      \
            FORT_FUNC_NAME(name, NAME, type_id, TYPE_ID)(FORT_WRAP_FUNC_ARG_2_VAL(lib, name));                 \
        }                                                                                                      \
    FORT_WRAP_EXTERN_C_END

#define FORT_WRAP_FUNC_DEF(lib, type, name, NAME)   \
    FORT_WRAP_FUNC_decl_and_def(lib, type, name, NAME, lib ## _ ## type ## _ ## type_id, lib ## _ ## type ## _ ## TYPE_ID)

#define FORT_WRAP_FUNC_DEF_2(lib, type, type_2, name, NAME)   \
    FORT_WRAP_FUNC_decl_and_def_2(lib, type, type_2, name, NAME, lib ## _ ## type ## _ ## type_id, lib ## _ ## type ## _ ## TYPE_ID)

#define FORT_WRAP_FUNC_DEF_real(lib, name, NAME) \
    FORT_WRAP_FUNC_DEF(lib, double, name, NAME)  \
    FORT_WRAP_FUNC_DEF(lib, float,  name, NAME)

#define FORT_WRAP_FUNC_DEF_complex(lib, name, NAME)      \
    FORT_WRAP_FUNC_DEF(lib, complex_double, name, NAME)  \
    FORT_WRAP_FUNC_DEF(lib, complex_float,  name, NAME)

#define FORT_WRAP_FUNC_DEF_all(lib, name, NAME)   \
    FORT_WRAP_FUNC_DEF_real(   lib, name, NAME)   \
    FORT_WRAP_FUNC_DEF_complex(lib, name, NAME)

#define FORT_WRAP_FUNC_DEF_2_complex(lib, name, NAME)              \
    FORT_WRAP_FUNC_DEF_2(lib, complex_double, double, name, NAME)  \
    FORT_WRAP_FUNC_DEF_2(lib, complex_float,  float,  name, NAME)

#endif /* FORT_FUNC_NAME */

/* -------------------------------------------------------------------------- */
