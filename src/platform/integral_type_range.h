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


#ifndef INTEGRAL_TYPE_RANGE_H
#define INTEGRAL_TYPE_RANGE_H

// -----------------------------------------------------------------------------

#include <climits>    // *_MIN and *_MAX values
#include <cstddef>    // std::{size_t, ptrdiff_t}

// -----------------------------------------------------------------------------
// Objective: Results of std::numeric_limits<>::min and max are known at
// compile time, but cannot be used for compile time checks because they are
// functions.  These classes make those limits accessible at compile time 
// (using the <climits> header) via a template argument.
// -----------------------------------------------------------------------------

template<typename T>
struct integral_type_range;

// -----------------------------------------------------------------------------

template<> struct integral_type_range<char>
{
    static const std::ptrdiff_t min_value = CHAR_MIN;
    static const std::size_t    max_value = CHAR_MAX;
    static const bool is_signed = true;
};

template<> struct integral_type_range<unsigned char>
{
    static const std::ptrdiff_t min_value = 0;
    static const std::size_t    max_value = UCHAR_MAX;
    static const bool is_signed = false;
};

template<> struct integral_type_range<short>
{
    static const std::ptrdiff_t min_value = SHRT_MIN;
    static const std::size_t    max_value = SHRT_MAX;
    static const bool is_signed = true;
};

template<> struct integral_type_range<unsigned short>
{
    static const std::ptrdiff_t min_value = 0;
    static const std::size_t    max_value = USHRT_MAX;
    static const bool is_signed = false;
};

template<> struct integral_type_range<int>
{
    static const std::ptrdiff_t min_value = INT_MIN;
    static const std::size_t    max_value = INT_MAX;
    static const bool is_signed = true;
};

template<> struct integral_type_range<unsigned int>
{
    static const std::ptrdiff_t min_value = 0;
    static const std::size_t    max_value = UINT_MAX;
    static const bool is_signed = false;
};

template<> struct integral_type_range<long>
{
    static const std::ptrdiff_t min_value = LONG_MIN;
    static const std::size_t    max_value = LONG_MAX;
    static const bool is_signed = true;
};

template<> struct integral_type_range<unsigned long>
{
    static const std::ptrdiff_t min_value = 0;
    static const std::size_t    max_value = ULONG_MAX;
    static const bool is_signed = false;
};

// We don't use std::ptrdiff_t and std::size_t types for
// min and max values in long long and unsigned long long
// specializations because they might not be sufficiently
// big to store them.

template<> struct integral_type_range<long long>
{
    static const long long          min_value = LLONG_MIN;
    static const unsigned long long max_value = LLONG_MAX;
    static const bool is_signed = true;
};

template<> struct integral_type_range<unsigned long long>
{
    static const long long          min_value = 0;
    static const unsigned long long max_value = ULLONG_MAX;
    static const bool is_signed = false;
};

// -----------------------------------------------------------------------------

// Check whether positive range of T1 is contained in positive range of T2.
template<typename T1, typename T2>
struct integral_type_range_check_positive
{
    static const bool value =
        integral_type_range<T1>::max_value <=
        integral_type_range<T2>::max_value;
};

// Check whether negative range of T1 is contained in negative range of T2.
template<typename T1, typename T2>
struct integral_type_range_check_negative
{
    static const bool value =
        integral_type_range<T2>::min_value <=
        integral_type_range<T1>::min_value;
};

// Check whether all T1 values are in the range of all T2 values.
template<typename T1, typename T2>
struct integral_type_range_check
{
    static const bool value =
        integral_type_range_check_positive<T1, T2>::value &&
        integral_type_range_check_negative<T1, T2>::value;
};

// -----------------------------------------------------------------------------

template<typename T2>
struct integral_type_range_check_val
{
    // Check whether t1 is in non-negative range of T2.
    template<typename T1>
    static bool in_non_negative_range(const T1& t1)
    {
        // Written in this way to avoid warnings due to using compile-time known
        // values in conditionals.

        if(t1 == 0)
        {
            return true;
        }
        else if(t1 > 0)
        {
            const std::size_t t11 = std::size_t(t1);
            return t11 <= integral_type_range<T2>::max_value;
        }
        else
        {
            return false;
        }
    }

    // Check whether t1 is in non-positive range of T2.
    template<typename T1>
    static bool in_non_positive_range(const T1& t1)
    {
        // Written in this way to avoid warnings due to using compile-time known
        // values in conditionals.

        if(t1 > 0)
        {
            return false;
        }
        else
        {
            const std::ptrdiff_t t11 = std::ptrdiff_t(t1);
            return t11 >= integral_type_range<T2>::min_value;
        }
    }

    // Check whether t1 is in range of T2.
    template<typename T1>
    static bool in_range(const T1& t1)
    {
        // Written in this way to avoid warnings due to using compile-time known
        // values in conditionals.

        if(t1 == 0)
        {
            return true;
        }
        else if(t1 > 0)
        {
            const std::size_t t11 = std::size_t(t1);
            return t11 <= integral_type_range<T2>::max_value;
        }
        else
        {
            const std::ptrdiff_t t11 = std::ptrdiff_t(t1);
            return t11 >= integral_type_range<T2>::min_value;
        }
    }
};

// -----------------------------------------------------------------------------

#endif // INTEGRAL_TYPE_RANGE_H

