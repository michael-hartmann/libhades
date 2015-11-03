#include <stdio.h>
#include <libhades.h>

int main(int argc, char *argv[])
{
    /* allocate a real 2x2 matrix */
    matrix_t *M = matrix_alloc(2,2);

    /* initialize matrix */
    matrix_set(M, 0,0, 1);
    matrix_set(M, 0,1, 2);
    matrix_set(M, 1,0, 3);
    matrix_set(M, 1,1, 4);

    printf("Inverse of\n");
    matrix_fprintf(stdout, M, "%+4g", "  ", "\n");

    /* invert matrix */
    matrix_invert(M);

    printf("is:\n");
    matrix_fprintf(stdout, M, "%+4g", "  ", "\n");

    /* free matrix */
    matrix_free(M);
     
    return 0;
}
