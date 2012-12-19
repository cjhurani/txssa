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


// test sparsity pattern for various inputs.

#include "p_norm_sparsity_vectors/p_norm_sparsity_vector.h"
#include "p_norm_of_vectors/p_norm_of_vectors.h"
#include "dense_vectors/dense_vectors.h"
#include <iostream>
#include <vector>

template<typename value_type>
const value_type& identity(const value_type& v)
{
    return v;
}

template<typename value_type, typename unary_func>
void vec_push_helper(
    value_type start,
    value_type stop,
    std::size_t n_steps,
    const unary_func& func,
    std::vector<value_type>& v)
{
    assert(start <= stop);

    const value_type step = (stop - start)/(n_steps - 1);

    for(std::size_t i = 0; i < n_steps; ++i)
    {
        v.push_back(func(start + i*step));
    }
}

template<typename value_type>
void test_p_norm_sparsity_vector(
    std::size_t& num_tests_done,
    std::size_t& num_tests_failed)
{
    std::vector<value_type> p_vals;
    std::vector<value_type> ratio_vals;
    std::vector<std::size_t> n_vals;

    p_vals.push_back(0);
    vec_push_helper<value_type, value_type (*)(value_type)>(-5.0, 0.0, 26, std::exp, p_vals);
    vec_push_helper<value_type, value_type (*)(value_type)>(0.1, 1.5, 20, std::exp, p_vals);
    p_vals.push_back(std::numeric_limits<value_type>::infinity());

    vec_push_helper<value_type>(0.0, 1.0, 201, identity<value_type>, ratio_vals);
    vec_push_helper<std::size_t>(1, 30, 30, identity<std::size_t>, n_vals);

    num_tests_done = 0;
    num_tests_failed = 0;

    std::vector<value_type> a_reserve(3*n_vals.back());
    std::vector<value_type> a_reserve_2(3*n_vals.back());
    std::vector<value_type> work(3*n_vals.back());
    std::vector<std::size_t> id(3*n_vals.back());

    for(std::size_t ip = 0; ip < p_vals.size(); ++ip)
    {
        const typename precision_traits<value_type>::scalar p = p_vals[ip];

        for(std::size_t ir = 0; ir < ratio_vals.size(); ++ir)
        {
            const value_type r = ratio_vals[ir];

            for(std::size_t in = 0; in < n_vals.size(); ++in)
            {
                const std::size_t n = n_vals[in];

                for(std::size_t j = 0; j < n; ++j)
                {
                    a_reserve[j] = ((value_type)rand())/RAND_MAX - 0.5;
                    if(a_reserve[j] == 0)
                        a_reserve[j] = 0.1;

                    a_reserve[j + n] = a_reserve[j];
                    a_reserve[j + 2*n] = a_reserve[j];
                }

                for(std::size_t m = 0; m < 3*n; ++m)
                {
                    ++num_tests_done;

                    std::size_t out_num_nnz;
                    bool success = p_norm_sparsity_vector<std::size_t, value_type>(
                        r,
                        p,
                        m,
                        3*n,
                        &a_reserve.front(),
                        1,
                        &out_num_nnz,
                        &id.front(),
                        &work.front());


                    if(num_tests_done % 1000 == 0)
                    {
                        std::cout << "Tests done: " << num_tests_done /1000 << "k " << " Tests failed: " << num_tests_failed << "\n";
                    }

                    if(success)
                    {
                        value_type p_norm, discarded_norm = 0;
                        dense_vectors<std::size_t, value_type> dv(1, 3*n, 3*n, &a_reserve.front());

                        p_norm_of_vectors(p, dv, &p_norm);

                        std::copy(a_reserve.begin(), a_reserve.begin() + 3*n, a_reserve_2.begin());

                        for(std::size_t ipat = 0; ipat < out_num_nnz; ++ipat)
                            a_reserve_2[id[ipat]] = 0;

                        dense_vectors<std::size_t, value_type> dv2(1, 3*n, 3*n, &a_reserve_2.front());
                        p_norm_of_vectors(p, dv2, &discarded_norm);

                        if(r == 1 && discarded_norm > 0)
                            success = false;

                        if(r < 1 && discarded_norm > (1 - r)*p_norm*(1 + 2*std::numeric_limits<value_type>::epsilon()))
                            success = false;

                        if(out_num_nnz < m)
                            success = false;
                    }
                    else
                    {
                        std::cout << "failed\n";
                        exit(1);
                    }

                    if(!success)
                        ++num_tests_failed;
                }
            }
        }
    }
}

int main()
{
    typedef double value_type;

    std::size_t num_tests_done, num_tests_failed;

    test_p_norm_sparsity_vector<value_type>(num_tests_done, num_tests_failed);

    std::cout << "num_tests_done   = " << num_tests_done << "\n";
    std::cout << "num_tests_failed = " << num_tests_failed << "\n";
}


