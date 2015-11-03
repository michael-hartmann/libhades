#include <stdio.h>
#include <libhades.h>
#include <libhades/expm.h>

int main(int argc, char *argv[])
{
    /* allocate memory for matrix M1 */
    matrix_complex_t *M = matrix_complex_alloc(3,3);

    /* initialize matrix:
     * This matrix is taken from
     * https://en.wikipedia.org/wiki/Matrix_exponential#Illustrations. This
     * matrix is not diagonalizable.
     */
    matrix_set(M, 0,0, 21);
    matrix_set(M, 0,1, 17);
    matrix_set(M, 0,2, 6);

    matrix_set(M, 1,0, -5);
    matrix_set(M, 1,1, -1);
    matrix_set(M, 1,2, -6);

    matrix_set(M, 2,0, 4);
    matrix_set(M, 2,1, 4);
    matrix_set(M, 2,2, 16);

    /* calculate matrix exponential */
    matrix_complex_expm(M);

    matrix_complex_fprintf(stdout, M, "%+5g%+0gi", "  ", "\n");

    /* free matrices M, B */
    matrix_complex_free(M);

    return 0;
}
