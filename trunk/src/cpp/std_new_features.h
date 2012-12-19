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


#ifndef STD_NEW_FEATURES_H
#define STD_NEW_FEATURES_H

// -----------------------------------------------------------------------------

// std_new_features_ prefix because these functions are already in C++11.
// Unfortunately, "g++ (GCC) 4.4.4 20100630 (Red Hat 4.4.4-10)" does not define
// __cplusplus to be 199711L when it should.  And there may be others.  So there
// is no guarantee that g++ will increment it for C++11.  So we duplicate the
// functions for all C++ standards, instead of using the C++11 std functions
// when using a C++11 compiler.

// -----------------------------------------------------------------------------

// Find extent that is sorted.  Returns a forward iterator set to the last
// element in sorted order.  For each value in sorted range, "next < current"
// is false.

template<typename val_it_fwd>
val_it_fwd std_new_features_is_sorted_until(
    val_it_fwd val_begin,
    val_it_fwd val_end)
{
    for(val_it_fwd next = val_begin;
        val_begin != val_end && ++next != val_end;
        ++val_begin)
    {
        if(*next < *val_begin)
            return val_begin;
    }

    return val_end;
}

// -----------------------------------------------------------------------------

// Test if range is sorted. If returns true, or each value in range,
// "next < current" is false.

template<typename val_it_fwd>
bool std_new_features_is_sorted(
    val_it_fwd val_begin,
    val_it_fwd val_end)
{
    return (std_new_features_is_sorted_until(val_begin, val_end) == val_end);
}

// -----------------------------------------------------------------------------

#endif // STD_NEW_FEATURES_H
