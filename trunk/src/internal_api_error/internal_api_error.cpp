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


// -----------------------------------------------------------------------------

#include "internal_api_error/internal_api_error.h"
#include <vector>
#include <string>
#include <limits>
#include <iostream>
#include <cstddef>

// -----------------------------------------------------------------------------

#ifdef _MSC_VER
#pragma warning( disable : 4514 ) // unreferenced inline function has been removed
#pragma warning( disable : 4710 ) // function not inlined
#endif

// -----------------------------------------------------------------------------

// anonymous
namespace
{
    std::vector<std::string> errors;
}

// -----------------------------------------------------------------------------

extern "C"
{

void internal_api_error_set_last(const char* str)
{
    try
    {
        errors.push_back(str);
    }
    catch(const std::exception& exc)
    {
        std::cerr
            << "internal_api_error_set_last: Exception. "
            << exc.what()
            << std::endl;
    }
    catch(...)
    {
        std::cerr
            << "internal_api_error_set_last: Exception. "
            << "Unknown"
            << std::endl;
    }
}

int internal_api_error_size()
{
    int ret = -1;
    const std::size_t size = errors.size();

    if(size <= std::size_t(std::numeric_limits<int>::max()))
        ret = int(errors.size());

    return ret;
}

int internal_api_error_string(int i, const char** ptr_to_error_string)
{
    int ret = -1;

    try
    {
        if(
            0 <= i &&
            std::size_t(i) < errors.size() &&
            ptr_to_error_string)
        {
            *ptr_to_error_string = errors[std::size_t(i)].c_str();
            ret = 0;
        }
    }
    catch(const std::exception& exc)
    {
        std::cerr
            << "internal_api_error_string: Exception. "
            << exc.what()
            << std::endl;
    }
    catch(...)
    {
        std::cerr
            << "internal_api_error_string: Exception. "
            << "Unknown"
            << std::endl;
    }

    return ret;
}

int internal_api_error_clear()
{
    int ret = -1;

    try
    {
        std::vector<std::string>().swap(errors);
        ret = 0;
    }
    catch(const std::exception& exc)
    {
        std::cerr
            << "internal_api_error_clear: Exception. "
            << exc.what()
            << std::endl;
    }
    catch(...)
    {
        std::cerr
            << "internal_api_error_clear: Exception. "
            << "Unknown"
            << std::endl;
    }

    return ret;
}

} // extern "C"

void internal_api_error_set_last(const std::string& str)
{
    internal_api_error_set_last(str.c_str());
}

// -----------------------------------------------------------------------------
