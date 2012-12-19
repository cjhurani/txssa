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


#include "txssa.h"

#include <string.h>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

// Objective: read a dense input matrix from a file and
// output a sparse spectral approximation based on other
// given arguments.
// Need in_file out_file value_type num_rows num_cols ...
// sparsity_ratio sparsity_norm_p max_num_bins

namespace {

template<typename T>
void read_val(std::ifstream& is, T& out)
{
    is >> out;
}

template<typename T>
void read_val(std::ifstream& is, std::complex<T>& out)
{
    T real, imag;
    is >> real >> imag;
    out = std::complex<T>(real, imag);
}

template<typename T>
void write_val(std::ofstream& os, const T& out)
{
    os << out;
}

template<typename T>
void write_val(std::ofstream& os, const std::complex<T>& out)
{
    os << out.real() << ' ' << out.imag();
}

} // namespace

template<typename value_type, typename scalar_type>
int ssa_file_io(
    const char* in_file,
    std::size_t num_rows,
    std::size_t num_cols,
    scalar_type sparsity_ratio,
    scalar_type sparsity_norm_p,
    std::size_t max_num_bins,
    const char* out_file)
{
    std::ifstream is(in_file);
    is.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    std::vector<value_type> col_values(num_rows*num_cols);

    for(std::size_t i = 0; i < num_rows; ++i)
        for(std::size_t j = 0; j < num_cols; ++j)
            read_val(is, col_values[i + j*num_rows]);

    bool impose_null_spaces = true;

    ssa_csr<std::size_t, std::size_t, value_type> out_matrix;

    int ret = ssa_lpn(
        num_rows, num_cols, &col_values.front(),
        num_rows, sparsity_ratio, sparsity_norm_p, max_num_bins,
        impose_null_spaces, ssa_matrix_type_general, out_matrix);

    if(ret != 0)
    {
        std::cerr << "ssa_lpn failed.\n";
        std::cerr << "Stack begin.\n";

        int n_errors = ssa_error_size();

        const char* ptr_to_error_string;

        for(int i = 0; i < n_errors; ++i)
        {
            if(ssa_error_string(i, &ptr_to_error_string) == 0)
                std::cerr << ptr_to_error_string << "\n";
            else
            {
                std::cerr << "Error in getting error stack\n";
            }
        }

        std::cerr << "Stack end.\n";

        ssa_error_clear();

        return 1;
    }

    std::ofstream os(out_file);
    os.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    os.precision(std::numeric_limits<scalar_type>::digits10);

    for(std::size_t i = 0; i < num_rows; ++i)
    {
        for(std::size_t jj = out_matrix.row_offsets[i]; jj < out_matrix.row_offsets[i+1]; ++jj)
        {
            os << i + 1 << ' ' << out_matrix.column_ids[jj] + 1 << ' ';
            write_val(os, out_matrix.values[jj]);
            os << '\n';
        }
    }

    return 0;
}


int main(int argc, char* argv[])
{
    if(argc < 9)
    {
        std::cerr << argv[0] << ": Insufficient number of arguments." << std::endl;
        std::cerr << "Need in_file out_file value_type num_rows num_cols sparsity_ratio sparsity_norm_p max_num_bins" << std::endl;
        exit(1);
    }

    const char* in_file      = argv[1];
    const char* out_file     = argv[2];
    const char* type_str     = argv[3];
    const char* num_rows_str = argv[4];
    const char* num_cols_str = argv[5];
    const char* ratio_str    = argv[6];
    const char* norm_p_str   = argv[7];
    const char* nbins_str    = argv[8];

    std::size_t num_rows = atoi(num_rows_str);
    std::size_t num_cols = atoi(num_cols_str);
    double ratio         = atof(ratio_str);
    double norm_p        = atof(norm_p_str);
    std::size_t nbins    = atoi(nbins_str);

    int ret = -1;

    try
    {
        if(!strcmp(type_str, "double"))
        {
            ret = ssa_file_io<double, double>(in_file, num_rows, num_cols,
                ratio, norm_p, nbins, out_file);
        }
        else if(!strcmp(type_str, "float"))
        {
            ret = ssa_file_io<float, float>(in_file, num_rows, num_cols,
                static_cast<float>(ratio), static_cast<float>(norm_p), nbins, out_file);
        }
        else if(!strcmp(type_str, "complex_double"))
        {
            ret = ssa_file_io< std::complex<double>, double>(in_file, num_rows, num_cols,
                ratio, norm_p, nbins, out_file);
        }
        else if(!strcmp(type_str, "complex_float"))
        {
            ret = ssa_file_io< std::complex<float>, float>(in_file, num_rows, num_cols,
                static_cast<float>(ratio), static_cast<float>(norm_p), nbins, out_file);
        }
        else
        {
            std::cerr << "Bad argument for type_str." << std::endl;
            exit(1);
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
        exit(1);
    }

    return ret;
}
