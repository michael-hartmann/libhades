#include <stdio.h>
#include <libhades.h>

int main(int argc, char *argv[])
{
    matrix_t *M = matrix_alloc(2,2);

    matrix_set(M, 0,0, 1);
    matrix_set(M, 0,1, 2);
    matrix_set(M, 1,0, 3);
    matrix_set(M, 1,1, 4);

    printf("Inverse of\n");
    matrix_fprintf(stdout, M, "%+4g", "  ", "\n");

    matrix_invert(M);

    printf("is:\n");
    matrix_fprintf(stdout, M, "%+4g", "  ", "\n");
     
    return 0;
}
