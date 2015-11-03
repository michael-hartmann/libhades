#include <stdio.h>
#include <libhades.h>

int main(int argc, char *argv[])
{
    matrix_t *M, *B;

    /* allocate memory for matrix M1 */
    M = matrix_alloc(3,3);

    /* initialize matrix */
    matrix_set(M, 0,0, 1);
    matrix_set(M, 0,1, 3);
    matrix_set(M, 0,2, -2);

    matrix_set(M, 1,0, 3);
    matrix_set(M, 1,1, 5);
    matrix_set(M, 1,2, 6);

    matrix_set(M, 2,0, 2);
    matrix_set(M, 2,1, 4);
    matrix_set(M, 2,2, 3);

    B = matrix_eye(3,NULL);

    matrix_solve(M, B);

    matrix_fprintf(stdout, B, "%+5g", "  ", "\n");

    /* free matrices M, B */
    matrix_free(M);
    matrix_free(B);

    return 0;
}
