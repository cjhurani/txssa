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

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

/*
 Objective: create a dense input matrix and compute two sparse spectral
 approximations
 1. from sparsity parameters
 2. using a given sparsity pattern.
*/

int main()
{
    const int n = 7;       /* matrix size */
    const int lda = n + 2; /* add 2 just for example (lda >= n) */

    const double sparsity_ratio = 0.8;
    const double sparsity_norm_p = 1.0;  /* L_p norm p value */
    const int max_num_bins = 256;        /* binning parameter */
    const int impose_null_spaces = 1;    /* impose if any exists */

    int i, j;
    int error;
    struct ssa_d_csr X1;

    /* allocate memory for dense matrix and the CSR pattern we want to give
       and the numerical matrix entries for that pattern */

    double* a               = (double*)malloc(sizeof(double)* (unsigned int)(n*lda));
    int* tri_diag_ids       = (int*)   malloc(sizeof(int)* (unsigned int)(3*n - 2));
    int* tri_diag_offsets   = (int*)   malloc(sizeof(int)* (unsigned int)(n + 1));
    double* tri_diag_values = (double*)malloc(sizeof(double)* (unsigned int)(3*n - 2));

    if(!a || !tri_diag_ids || !tri_diag_offsets || !tri_diag_values)
    {
        fprintf(stderr, "%s\n", "malloc failed.");
        exit(1);
    }

    /* create a (non-strict) diagonally dominant symmetric singular matrix */

    for(j = 0; j < n; ++j)
        for(i = 0; i < n; ++i)
            a[i + lda*j] = i == j ? 0 : - (i + j);

    for(j = 0; j < n; ++j)
        for(i = 0; i < n; ++i)
            if(i != j)
                a[j + lda*j] -= a[i + lda*j];

    /* compute X1 (the sparse approximation) without giving its pattern.
       one has to give sparsity_ratio and sparsity_norm_p */

    error = ssa_d_lpn(
        n, n, a, lda, sparsity_ratio, sparsity_norm_p,
        max_num_bins, impose_null_spaces, ssa_matrix_type_hermitian,
        &X1);

    if(error)
    {
        fprintf(stderr, "%s\n", "ssa_d_lpn failed.");
        exit(1);
    }

    /* Now we create and use a tridiagonal input pattern */

    tri_diag_offsets[0] = 0;

    for(i = 0; i < n; ++i)
        tri_diag_offsets[i + 1] = tri_diag_offsets[i] + ((i == 0 || i == n - 1) ? 2 : 3);

    if(n == 1)
    {
        tri_diag_ids[0] = 0;
    }
    else if(n == 2)
    {
        tri_diag_ids[0] = 0;
        tri_diag_ids[1] = 1;
        tri_diag_ids[2] = 0;
        tri_diag_ids[3] = 1;
    }
    else
    {
        tri_diag_ids[0] = 0;
        tri_diag_ids[1] = 1;

        for(i = 1; i < n - 1; ++i)
        {
            tri_diag_ids[2 + 3*(i-1)] = i-1;
            tri_diag_ids[3 + 3*(i-1)] = i;
            tri_diag_ids[4 + 3*(i-1)] = i+1;
        }

        tri_diag_ids[3*(n-1) - 1] = n-2;
        tri_diag_ids[3*(n-1)] = n-1;
    }

    /* compute the sparse approximation values using the tridiagonal pattern.
       one does not have to give sparsity_ratio and sparsity_norm_p */

    error = ssa_d_pat(
        n, n, a, lda, tri_diag_offsets, tri_diag_ids,
        max_num_bins, impose_null_spaces, ssa_matrix_type_hermitian,
        tri_diag_values);

    if(error)
    {
        fprintf(stderr, "%s\n", "ssa_d_pat failed.");
        exit(1);
    }

    printf("%s = %d\n", "ssa_impl_version", ssa_impl_version);

    printf("A:\n");
    for(i = 0; i < n; ++i)
    {
        for(j = 0; j < n; ++j)
        {
            printf("%10.4f ", a[i + lda*j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("%s\n", "X1 using");
    printf("%s %7.4f\n", "sparsity_ratio  =   ", sparsity_ratio);
    printf("%s %7.4f\n", "sparsity_norm_p =   ", sparsity_norm_p);
    printf("%s %d\n",    "max_num_bins =       ", max_num_bins);
    printf("%s %d\n",    "impose_null_spaces = ", impose_null_spaces);

    for(i = 0; i < n; ++i)
    {
        printf("%s %d: ", "row = ", i);
        for(j = X1.row_offsets[i]; j < X1.row_offsets[i + 1]; ++j)
            printf("(%d, %g) ", X1.column_ids[j], X1.values[j]);

        printf("\n");
    }
    printf("\n");

    printf("%s\n", "X2 using tridiagonal input pattern and");
    printf("%s %d\n", "max_num_bins =       ", max_num_bins);
    printf("%s %d\n", "impose_null_spaces = ", impose_null_spaces);

    for(i = 0; i < n; ++i)
    {
        printf("%s %d: ", "row = ", i);
        for(j = tri_diag_offsets[i]; j < tri_diag_offsets[i + 1]; ++j)
            printf("(%d, %g) ", tri_diag_ids[j], tri_diag_values[j]);

        printf("\n");
    }

    ssa_d_csr_deallocate(&X1);
    free(tri_diag_values);
    free(tri_diag_offsets);
    free(tri_diag_ids);
    free(a);

    return 0;
}
