#include <math.h>
#include <stdio.h>
#include <complex.h>

#include <lapacke.h>
#include <libhades.h>

#include "unittest.h"

int test_trace()
{
    matrix_t *A, *B;
    matrix_complex_t *C, *D;
    double trace;
    complex_t trace_cplx;
    unittest_t test;
    unittest_init(&test, "matrix_mult", "Test trace");

    A = matrix_alloc(3,3);
    B = matrix_alloc(3,3);

    matrix_set(A, 0,0, 1);
    matrix_set(A, 0,1, 2);
    matrix_set(A, 0,2, 3);
    matrix_set(A, 1,0, 4);
    matrix_set(A, 1,1, 5);
    matrix_set(A, 1,2, 6);
    matrix_set(A, 2,0, 7);
    matrix_set(A, 2,1, 8);
    matrix_set(A, 2,2, 9);

    matrix_set(B, 0,0, 9);
    matrix_set(B, 0,1, 8);
    matrix_set(B, 0,2, 7);
    matrix_set(B, 1,0, 6);
    matrix_set(B, 1,1, 5);
    matrix_set(B, 1,2, 4);
    matrix_set(B, 2,0, 3);
    matrix_set(B, 2,1, 2);
    matrix_set(B, 2,2, 1);

    trace = matrix_trace_AB(A, B);
    AssertEqual(&test, trace, 189);

    AssertEqual(&test, matrix_trace(A), 15);
    AssertEqual(&test, matrix_trace(B), 15);

    matrix_free(A);
    matrix_free(B);

    C = matrix_complex_alloc(3,3);
    D = matrix_complex_alloc(3,3);

    matrix_set(C, 0,0, CPLX(0,1));
    matrix_set(C, 0,1, CPLX(-2,0));
    matrix_set(C, 0,2, CPLX(3,1));
    matrix_set(C, 1,0, CPLX(-2,0));
    matrix_set(C, 1,1, CPLX(1,0));
    matrix_set(C, 1,2, CPLX(0,1));
    matrix_set(C, 2,0, CPLX(0,-2));
    matrix_set(C, 2,1, CPLX(1,3));
    matrix_set(C, 2,2, CPLX(5,0));

    matrix_set(D, 0,0, CPLX(0,-2));
    matrix_set(D, 0,1, CPLX(0,1));
    matrix_set(D, 0,2, CPLX(2,-2));
    matrix_set(D, 1,0, CPLX(1,0));
    matrix_set(D, 1,1, CPLX(5,0));
    matrix_set(D, 1,2, CPLX(0,1));
    matrix_set(D, 2,0, CPLX(-1,0));
    matrix_set(D, 2,1, CPLX(-4,-4));
    matrix_set(D, 2,2, CPLX(0,1));

    trace_cplx = matrix_complex_trace(C);
    AssertEqual(&test, CREAL(trace_cplx), 6);
    AssertEqual(&test, CIMAG(trace_cplx), 1);
    trace_cplx = matrix_complex_trace(D);
    AssertEqual(&test, CREAL(trace_cplx), 5);
    AssertEqual(&test, CIMAG(trace_cplx), -1);

    trace = matrix_trace_complex_AB_real(C, D);
    AssertEqual(&test, trace, -1);

    trace_cplx = matrix_trace_complex_AB(C, D);
    AssertEqual(&test, CREAL(trace_cplx), -1);
    AssertEqual(&test, CIMAG(trace_cplx), -5);

    matrix_complex_free(C);
    matrix_complex_free(D);

    return test_results(&test, stderr);
}

int test_matrix_mult()
{
    matrix_t *A, *B, *C;
    unittest_t test;
    unittest_init(&test, "matrix_mult", "Test matrix multiplication");

    A = matrix_alloc(2,3);
    B = matrix_alloc(3,2);

    matrix_set(A, 0,0, 3);
    matrix_set(A, 0,1, 2);
    matrix_set(A, 0,2, 1);
    matrix_set(A, 1,0, 1);
    matrix_set(A, 1,1, 0);
    matrix_set(A, 1,2, 2);

    matrix_set(B, 0,0, 1);
    matrix_set(B, 0,1, 2);
    matrix_set(B, 1,0, 0);
    matrix_set(B, 1,1, 1);
    matrix_set(B, 2,0, 4);
    matrix_set(B, 2,1, 0);

    C = matrix_mult(A,B,1,NULL);

    matrix_free(A);
    matrix_free(B);

    AssertEqual(&test, matrix_get(C,0,0), 7);
    AssertEqual(&test, matrix_get(C,0,1), 8);
    AssertEqual(&test, matrix_get(C,1,0), 9);
    AssertEqual(&test, matrix_get(C,1,1), 2);

    matrix_free(C);

    return test_results(&test, stderr);
}

int test_matrix_zeros()
{
    int m,n;
    matrix_t *A;
    matrix_complex_t *C;
    unittest_t test;

    unittest_init(&test, "matrix_zeros", "Test function matrix_zeros");

    for(m = 6; m < 20; m++)
        for(n = 1; n < 8; n++)
        {
            int im, in;
            A = matrix_zeros(m,n, NULL);
            C = matrix_complex_zeros(m,n, NULL);
            for(im = 0; im < m; im++)
                for(in = 0; in < n; in++)
                {
                    complex_t c = matrix_get(C,im,in);
                    AssertEqual(&test, matrix_get(A,im,in),0);
                    AssertEqual(&test, CREAL(c), 0);
                    AssertEqual(&test, CIMAG(c), 0);
                }

            matrix_free(A);
            matrix_complex_free(C);
        }


    return test_results(&test, stderr);
}

int test_zeros()
{
    matrix_t *A;
    matrix_complex_t *B;
    unittest_t test;
    unittest_init(&test, "matrix_zeros", "Test zeros");
    size_t rows = 100, columns = 150;

    A = matrix_zeros        (rows, columns, NULL);
    B = matrix_complex_zeros(rows, columns, NULL);

    for(size_t i = 0; i < rows; i++)
        for(size_t j = 0; j < columns; j++)
        {
            AssertEqual(&test, matrix_get(A, i,j), 0);
            AssertEqual(&test, matrix_get(B, i,j), 0);
        }

    matrix_free(A);
    matrix_complex_free(B);

    return test_results(&test, stderr);
}

int test_matrix_kron()
{
    matrix_t *A, *B, *AB;
    matrix_complex_t *C, *D, *CD;
    matrix_t *E, *EE;
    unittest_t test;
    unittest_init(&test, "matrix_kron", "Test function matrix_kron");

    A = matrix_alloc(3,2);
    B = matrix_alloc(2,2);

    matrix_set(A, 0,0, 1);
    matrix_set(A, 0,1, 2);
    matrix_set(A, 1,0, 3);
    matrix_set(A, 1,1, 4);
    matrix_set(A, 2,0, 5);
    matrix_set(A, 2,1, 6);

    matrix_set(B, 0,0, 7);
    matrix_set(B, 0,1, 8);
    matrix_set(B, 1,0, 9);
    matrix_set(B, 1,1, 0);

    AB = matrix_kron(A, B, NULL);

    matrix_free(A);
    matrix_free(B);

    AssertEqual(&test, matrix_get(AB, 0,0), 7);
    AssertEqual(&test, matrix_get(AB, 0,1), 8);
    AssertEqual(&test, matrix_get(AB, 0,2), 14);
    AssertEqual(&test, matrix_get(AB, 0,3), 16);

    AssertEqual(&test, matrix_get(AB, 1,0), 9);
    AssertEqual(&test, matrix_get(AB, 1,1), 0);
    AssertEqual(&test, matrix_get(AB, 1,2), 18);
    AssertEqual(&test, matrix_get(AB, 1,3), 0);

    AssertEqual(&test, matrix_get(AB, 2,0), 21);
    AssertEqual(&test, matrix_get(AB, 2,1), 24);
    AssertEqual(&test, matrix_get(AB, 2,2), 28);
    AssertEqual(&test, matrix_get(AB, 2,3), 32);

    AssertEqual(&test, matrix_get(AB, 3,0), 27);
    AssertEqual(&test, matrix_get(AB, 3,1), 0);
    AssertEqual(&test, matrix_get(AB, 3,2), 36);
    AssertEqual(&test, matrix_get(AB, 3,3), 0);

    AssertEqual(&test, matrix_get(AB, 4,0), 35);
    AssertEqual(&test, matrix_get(AB, 4,1), 40);
    AssertEqual(&test, matrix_get(AB, 4,2), 42);
    AssertEqual(&test, matrix_get(AB, 4,3), 48);

    AssertEqual(&test, matrix_get(AB, 5,0), 45);
    AssertEqual(&test, matrix_get(AB, 5,1), 0);
    AssertEqual(&test, matrix_get(AB, 5,2), 54);
    AssertEqual(&test, matrix_get(AB, 5,3), 0);

    matrix_free(AB);

    C = matrix_complex_alloc(2,2);
    D = matrix_complex_alloc(2,2);

    matrix_set(C, 0,0, CPLX(1,1));
    matrix_set(C, 0,1, CPLX(2,-2));
    matrix_set(C, 1,0, CPLX(5,0));
    matrix_set(C, 1,1, CPLX(0,-1));

    matrix_set(D, 0,0, CPLX(0,-1));
    matrix_set(D, 0,1, CPLX(2,0));
    matrix_set(D, 1,0, CPLX(5,0));
    matrix_set(D, 1,1, CPLX(0,7));

    CD = matrix_complex_alloc(4,4);
    matrix_complex_kron(C, D, CD);

    matrix_complex_free(C);
    matrix_complex_free(D);

    AssertEqual(&test, CREAL(matrix_get(CD, 0,0)), 1);
    AssertEqual(&test, CIMAG(matrix_get(CD, 0,0)), -1);
    AssertEqual(&test, CREAL(matrix_get(CD, 0,1)), 2);
    AssertEqual(&test, CIMAG(matrix_get(CD, 0,1)), 2);
    AssertEqual(&test, CREAL(matrix_get(CD, 0,2)), -2);
    AssertEqual(&test, CIMAG(matrix_get(CD, 0,2)), -2);
    AssertEqual(&test, CREAL(matrix_get(CD, 0,3)), 4);
    AssertEqual(&test, CIMAG(matrix_get(CD, 0,3)), -4);

    AssertEqual(&test, CREAL(matrix_get(CD, 1,0)), 5);
    AssertEqual(&test, CIMAG(matrix_get(CD, 1,0)), 5);
    AssertEqual(&test, CREAL(matrix_get(CD, 1,1)), -7);
    AssertEqual(&test, CIMAG(matrix_get(CD, 1,1)), 7);
    AssertEqual(&test, CREAL(matrix_get(CD, 1,2)), 10);
    AssertEqual(&test, CIMAG(matrix_get(CD, 1,2)), -10);
    AssertEqual(&test, CREAL(matrix_get(CD, 1,3)), 14);
    AssertEqual(&test, CIMAG(matrix_get(CD, 1,3)), 14);

    AssertEqual(&test, CREAL(matrix_get(CD, 2,0)), 0);
    AssertEqual(&test, CIMAG(matrix_get(CD, 2,0)), -5);
    AssertEqual(&test, CREAL(matrix_get(CD, 2,1)), 10);
    AssertEqual(&test, CIMAG(matrix_get(CD, 2,1)), 0);
    AssertEqual(&test, CREAL(matrix_get(CD, 2,2)), -1);
    AssertEqual(&test, CIMAG(matrix_get(CD, 2,2)), 0);
    AssertEqual(&test, CREAL(matrix_get(CD, 2,3)), 0);
    AssertEqual(&test, CIMAG(matrix_get(CD, 2,3)), -2);

    AssertEqual(&test, CREAL(matrix_get(CD, 3,0)), 25);
    AssertEqual(&test, CIMAG(matrix_get(CD, 3,0)), 0);
    AssertEqual(&test, CREAL(matrix_get(CD, 3,1)), 0);
    AssertEqual(&test, CIMAG(matrix_get(CD, 3,1)), 35);
    AssertEqual(&test, CREAL(matrix_get(CD, 3,2)), 0);
    AssertEqual(&test, CIMAG(matrix_get(CD, 3,2)), -5);
    AssertEqual(&test, CREAL(matrix_get(CD, 3,3)), 7);
    AssertEqual(&test, CIMAG(matrix_get(CD, 3,3)), 0);

    matrix_complex_free(CD);

    E = matrix_zeros(2,2, NULL);

    matrix_set(E, 0,0, 1);
    matrix_set(E, 1,1, 2);

    EE = matrix_kron(E, E, NULL);
    matrix_free(E);

    AssertEqual(&test, matrix_get(EE, 0,0), 1);
    AssertEqual(&test, matrix_get(EE, 1,1), 2);
    AssertEqual(&test, matrix_get(EE, 2,2), 2);
    AssertEqual(&test, matrix_get(EE, 3,3), 4);

    for(size_t i = 0; i < 4; i++)
        for(size_t j = 0; j < 4; j++)
            if(i != j)
                AssertEqual(&test, matrix_get(EE, i,j), 0);


    matrix_free(EE);

    return test_results(&test, stderr);
}


int main(int argc, char *argv[])
{
    test_matrix_zeros();
    test_matrix_kron();
    test_matrix_mult();
    test_zeros();
    test_trace();
    
    return 0;
}
