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


#ifndef COMPLEX_TYPES_H
#define COMPLEX_TYPES_H

/* -------------------------------------------------------------------------- */

/* This header can be included in a C translation unit. */

#ifdef __cplusplus

#include <complex>

typedef std::complex<float>  complex_float;
typedef std::complex<double> complex_double;

#else

typedef struct complex_float_tag
{
    float real;
    float imag;
} complex_float;

typedef struct complex_double_tag
{
    double real;
    double imag;
} complex_double;

#endif

/* -------------------------------------------------------------------------- */

#ifdef __cplusplus

// -----------------------------------------------------------------------------

namespace std
{
    // Functions take const& to match the complex versions.

    inline float real(const float& a)
    {
        return a;
    }

    inline double real(const double& a)
    {
        return a;
    }

    inline float imag(const float& /*a*/)
    {
        return 0;
    }

    inline double imag(const double& /*a*/)
    {
        return 0;
    }

    inline float conj(const float& a)
    {
        return a;
    }

    inline double conj(const double& a)
    {
        return a;
    }

    inline double abs_square(const double& a)
    {
        return a*a;
    }

    inline float abs_square(const float& a)
    {
        return a*a;
    }

    inline double abs_square(const complex_double& a)
    {
        return a.real() * a.real() + a.imag() * a.imag();
    }

    inline float abs_square(const complex_float& a)
    {
        return a.real() * a.real() + a.imag() * a.imag();
    }
}

// -----------------------------------------------------------------------------

// complex_rebind<new_type, float>::other == new_type
// but
// complex_rebind<new_type, std::complex<float> >::other == std::complex<new_type>
//
// And so on on for scalar == double.

template<typename new_type, typename scalar>
struct complex_rebind
{
    typedef new_type other;
};

template<typename new_type, typename scalar>
struct complex_rebind<new_type, std::complex<scalar> >
{
    typedef std::complex<new_type> other;
};

// -----------------------------------------------------------------------------

// Extractor functors.  Inlinable compared to std::{real,imag} for std::complex.

template<typename scalar_type>
struct real_extractor
{
    inline scalar_type operator()(const std::complex<scalar_type>& a) const
    {
        return a.real();
    }

    inline scalar_type operator()(scalar_type a) const
    {
        return a;
    }
};

template<typename scalar_type>
struct imag_extractor
{
    inline scalar_type operator()(const std::complex<scalar_type>& a) const
    {
        return a.imag();
    }

    inline scalar_type operator()(scalar_type a) const
    {
        return 0;
    }
};

template<typename scalar_type>
struct conj_imag_extractor
{
    inline scalar_type operator()(const std::complex<scalar_type>& a) const
    {
        return - a.imag();
    }

    inline scalar_type operator()(scalar_type a) const
    {
        return 0;
    }
};

// -----------------------------------------------------------------------------

#endif /* __cplusplus */

/* -------------------------------------------------------------------------- */

#endif /* COMPLEX_TYPES_H */
