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


#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

// -----------------------------------------------------------------------------

#include <cstddef>
#include <cmath>       // std::abs
#include <cassert>

// -----------------------------------------------------------------------------
// Objective: Small utility functions for working on vectors.
// -----------------------------------------------------------------------------

// Compute inf norm of difference of vector and a value.  Ignore those entries
// in the vector that are equal to "to_ignore".  It is assumed that floating
// point issues don't happen in the ignoring check.
template<typename index_type, typename value_type>
value_type vector_utils_max_abs_diff_with_ignore(
    index_type size,
    const value_type* vec,
    value_type val,
    value_type to_ignore)
{
    assert(vec);

    value_type max_abs_diff = 0;

    for(index_type i = 0; i < size; ++i)
    {
        if(vec[i] != to_ignore)
        {
            const value_type abs_diff = std::abs(vec[i] - val);

            if( max_abs_diff < abs_diff)
                max_abs_diff = abs_diff;
        }
    }

    return max_abs_diff;
}

// -----------------------------------------------------------------------------

// Invert non-zero values.
template<typename index_type, typename value_type>
void vector_utils_invert_non_zero(
    index_type size,
    value_type* vec)
{
    assert(vec);

    for(index_type i = 0; i < size; ++i)
    {
        if(vec[i] != 0)
        {
            vec[i] = 1/vec[i];
        }
    }
}

// -----------------------------------------------------------------------------

// Replace each value with its p-th power.
template<typename index_type, typename value_type, typename pow_type>
void vector_utils_replace_with_pow(
    index_type size,
    value_type* vec,
    const pow_type& p)
{
    assert(vec);

    for(index_type i = 0; i < size; ++i)
    {
        vec[i] = std::pow(vec[i], p);
    }
}

// -----------------------------------------------------------------------------

// y <- a * x + b * y

template<typename index_type, typename value_type>
void vector_utils_axpby(
    index_type n,
    const value_type* x, // can be 0 if a is value_type(0)
    value_type* y,
    value_type a,
    value_type b,
    index_type incx,
    index_type incy)
{
    assert(a == value_type(0) || x);
    assert(y);
    assert(0 < incx);
    assert(0 < incy);

    // y <- a * x + b * y

    // a = -1, b = -1 => y = - x - y
    // a = -1, b =  0 => y =  -x
    // a = -1, b =  1 => y =  y - x

    // a =  0, b = -1 => y = -y
    // a =  0, b =  0 => y =  0
    // a =  0, b =  1 => y =  y

    // a =  1, b = -1 => y = x - y
    // a =  1, b =  0 => y = x
    // a =  1, b =  1 => y = x + y

    // a = -1         => y = b * y - x
    // a =  0         => y = b * y
    // a =  1         => y = b * y + x

    //         b = -1 => y = a * x - y
    //         b =  0 => y = a * x
    //         b =  1 => y = a * x + y

    const std::size_t nx = std::size_t(n) * std::size_t(incx);

    if(a == value_type(-1))
    {
        if(b == value_type(-1))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = - y[iy] - x[ix];
            return;
        }
        if(b == value_type(0))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = - x[ix];
            return;
        }
        if(b == value_type(1))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] -= x[ix];
            return;
        }

        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] = b * y[iy] - x[ix];
        return;
    }

    if(a == value_type(0))
    {
        if(b == value_type(-1))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = -y[iy];
            return;
        }
        if(b == value_type(0))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = 0;
            return;
        }
        if(b == value_type(1))
        {
            // y = y
            return;
        }

        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] *= b;
        return;
    }

    if(a == value_type(1))
    {
        if(b == value_type(-1))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = x[ix] - y[iy];
            return;
        }
        if(b == value_type(0))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] = x[ix];
            return;
        }
        if(b == value_type(1))
        {
            for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
                y[iy] += x[ix];
            return;
        }

        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] = b * y[iy] + x[ix];
        return;
    }

    if(b == value_type(-1))
    {
        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] = a*x[ix] - y[iy];
        return;
    }
    if(b == value_type(0))
    {
        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] = a*x[ix];
        return;
    }
    if(b == value_type(1))
    {
        for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
            y[iy] += a*x[ix];
        return;
    }

    for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
        y[iy] = a*x[ix] + b*y[iy];
}

// -----------------------------------------------------------------------------

template<typename index_type, typename value_type>
void vector_utils_add(
    index_type n,
    const value_type* x,
    value_type* y,
    index_type incx,
    index_type incy)
{
    assert(x);
    assert(y);
    assert(0 < incx);
    assert(0 < incy);

    const std::size_t nx = std::size_t(n) * std::size_t(incx);

    for(std::size_t ix = 0, iy = 0; ix < nx; ix += incx, iy += incy)
        y[iy] += x[ix];
}

// -----------------------------------------------------------------------------

#endif // VECTOR_UTILS_H
