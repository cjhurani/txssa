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


#ifndef CPU_TIMER_H
#define CPU_TIMER_H

// -----------------------------------------------------------------------------

#ifdef _MSC_VER
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <ctime>
#endif

#include <cassert>
#include <limits>

// -----------------------------------------------------------------------------

class cpu_timer
{
    // -------------------------------------------------------------------------

public:

    cpu_timer()
        : delta_t(std::numeric_limits<double>::quiet_NaN()),
          finished(false),
          started(false)
    {
#ifdef _MSC_VER
        QueryPerformanceCounter(&start);
#else
        timespec ts;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
        start = ts.tv_sec + 1E-9 * ts.tv_nsec;
#endif
        started = true;
    }

    ~cpu_timer()
    {
        done();
    }

    double done()
    {
        if(started && !finished)
        {
#ifdef _MSC_VER
            QueryPerformanceCounter(&stop);

            LARGE_INTEGER freq;
            QueryPerformanceFrequency(&freq);

            delta_t = double(stop.QuadPart - start.QuadPart) / freq.QuadPart;
#else
            timespec ts;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts);
            stop = ts.tv_sec + 1E-9 * ts.tv_nsec;
            delta_t = stop - start;
#endif
            finished = true;
        }

        return delta_t;
    }

    double time_taken() const
    {
        assert((started && !finished) == false);

        return delta_t;
    }

private:

#ifdef _MSC_VER
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
#else
    double start;
    double stop;
#endif

    double delta_t;
    bool finished;
    bool started;
};

// -----------------------------------------------------------------------------

#endif // CPU_TIMER_H
