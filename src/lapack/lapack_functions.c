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


#include "lapack/lapack_functions.h"
#include "fortran/fort_wrap_func_def.h"

FORT_WRAP_FUNC_DEF_all(      LAPACK, trtri, TRTRI)
FORT_WRAP_FUNC_DEF_all(      LAPACK, geqrf, GEQRF)
FORT_WRAP_FUNC_DEF_all(      LAPACK, potrf, POTRF)
FORT_WRAP_FUNC_DEF_all(      LAPACK, posv,  POSV)
FORT_WRAP_FUNC_DEF_real(     LAPACK, orgqr, ORGQR)
FORT_WRAP_FUNC_DEF_complex(  LAPACK, ungqr, UNGQR)
FORT_WRAP_FUNC_DEF_real(     LAPACK, ormqr, ORMQR)
FORT_WRAP_FUNC_DEF_complex(  LAPACK, unmqr, UNMQR)
FORT_WRAP_FUNC_DEF_real(     LAPACK, geqp3, GEQP3)
FORT_WRAP_FUNC_DEF_2_complex(LAPACK, geqp3, GEQP3)
FORT_WRAP_FUNC_DEF_real(     LAPACK, gesvd, GESVD)
FORT_WRAP_FUNC_DEF_2_complex(LAPACK, gesvd, GESVD)
