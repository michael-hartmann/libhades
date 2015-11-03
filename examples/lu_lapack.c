#include <stdio.h>
#include <libhades.h>
#include <libhades/expm.h>

/* This example illustrates how to access LAPACK functions directly. */
int main(int argc, char *argv[])
{
    int info = -1;
    int ipiv[3];
    /* allocate memory for matrix M */
    matrix_t *M = matrix_alloc(3,3);

    /* initialize matrix: */
    matrix_set(M, 0,0, 21);
    matrix_set(M, 0,1, 17);
    matrix_set(M, 0,2, 6);

    matrix_set(M, 1,0, -5);
    matrix_set(M, 1,1, -1);
    matrix_set(M, 1,2, -6);

    matrix_set(M, 2,0, 4);
    matrix_set(M, 2,1, 4);
    matrix_set(M, 2,2, 16);

    /* Note:
     * 1) LAPACK function names have a underscore at the end
     * 2) We could also use the libhades function matrix_lu_decomposition(M)
     *    here
     */
    dgetrf_(
        &M->rows,    /* M number of rows of M */
        &M->columns, /* N number of columns of M */
        M->M,        /* matrix M to be factored */
        &M->columns, /* LDA: leading dimension of M */
        ipiv,        /* pivot indices of dimension (min(M,N)) */
        &info
    );

    printf("info = %d\n", info);
    matrix_fprintf(stdout, M, "%g", "  ", "\n");

    /* free matrices M, B */
    matrix_free(M);

    return 0;
}
