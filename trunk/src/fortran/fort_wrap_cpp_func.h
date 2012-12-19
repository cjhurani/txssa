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


#ifndef FORT_WRAP_CPP_FUNC_H
#define FORT_WRAP_CPP_FUNC_H

/* This header can be included in a C translation unit. */

/* This file should be included by the caller of the wrapper functions. */

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus

// -----------------------------------------------------------------------------

#include "fortran/fort_wrap_func.h"

// -----------------------------------------------------------------------------

// inline for weak symbols.

#define FORT_WRAP_CPP_FUNC_DEF(lib, type, name)                              \
    inline FORT_RET lib ## _ ## name(lib ## _ ## name ## _ARG(type))         \
    {                                                                        \
        FORT_WRAP_FUNC(lib, type, name)(lib ## _ ## name ## _ARG_VAL);       \
    }

#define FORT_WRAP_CPP_FUNC_DEF_ARG_2(lib, type_1, type_2, name)              \
    inline                                                                   \
    FORT_RET lib ## _ ## name(lib ## _ ## name ## _ARG_2(type_1, type_2))    \
    {                                                                        \
        FORT_WRAP_FUNC(lib, type_1, name)(lib ## _ ## name ## _ARG_2_VAL);   \
    }

#define FORT_WRAP_CPP_FUNC_DEF_real(lib, name)   \
    FORT_WRAP_CPP_FUNC_DEF(lib, float,  name)    \
    FORT_WRAP_CPP_FUNC_DEF(lib, double, name)

#define FORT_WRAP_CPP_FUNC_DEF_complex(lib, name)      \
    FORT_WRAP_CPP_FUNC_DEF(lib, complex_float, name)   \
    FORT_WRAP_CPP_FUNC_DEF(lib, complex_double, name)

#define FORT_WRAP_CPP_FUNC_DEF_all(lib, name)   \
    FORT_WRAP_CPP_FUNC_DEF_complex(lib, name)   \
    FORT_WRAP_CPP_FUNC_DEF_real(lib, name)

#define FORT_WRAP_CPP_FUNC_DEF_ARG_2_complex(lib, name)              \
    FORT_WRAP_CPP_FUNC_DEF_ARG_2(lib, complex_float,  float, name)   \
    FORT_WRAP_CPP_FUNC_DEF_ARG_2(lib, complex_double, double, name)

// -----------------------------------------------------------------------------

#endif /* __cplusplus */

/* -------------------------------------------------------------------------- */

#endif /* FORT_WRAP_CPP_FUNC_H */
