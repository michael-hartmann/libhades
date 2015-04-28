/**
 * @file   expm.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   April, 2015
 * @brief  implement matrix exponential function
 */


#include <math.h>
#include <libhades.h>
#include <libhades/expm.h>

#define THETA_3  1.5e-2
#define THETA_5  2.5e-1
#define THETA_7  9.5e-1
#define THETA_9  2.1e0
#define THETA_13 5.4e0

static double const b3[]  = {120,60,12,1};
static double const b5[]  = {30240,15120,3360,420,30,1};
static double const b7[]  = {17297280, 8648640,1995840, 277200, 25200,1512, 56,1};
static double const b9[]  = {17643225600,8821612800,2075673600,302702400,30270240, 2162160,110880,3960,90,1};
static double const b13[] = {64764752532480000.,32382376266240000.,7771770303897600.,1187353796428800.,129060195264000.,10559470521600.,670442572800.,33522128640.,1323241920.,40840800.,960960.,16380.,182,1};

static double const *b_list[] = {
    NULL, /* 0 */
    NULL, /* 1 */
    NULL, /* 2 */
    b3,   /* 3 */
    NULL, /* 4 */
    b5,   /* 5 */
    NULL, /* 6 */
    b7,   /* 7 */
    NULL, /* 8 */
    b9,   /* 9 */
    NULL, /* 10 */
    NULL, /* 11 */
    NULL, /* 12 */
    b13   /* 13 */
};

/* works only for square matrices and M == 3,5,7 */
static matrix_complex_t *_expm_pade3579(matrix_complex_t *A, int M)
{
    matrix_complex_t *U = NULL, *V = NULL, *A2 = NULL, *A2n = NULL, *X;
    const int dim = A->rows;
    double const *b = b_list[M];

    /* set U = b[1]*Id and V = b[0]*Id */
    U = matrix_complex_eye(dim, NULL);
    V = matrix_complex_eye(dim, NULL);
    for(int i = 0; i < dim; i++)
    {
        matrix_set(U, i,i, matrix_get(U,i,i)*b[1]);
        matrix_set(V, i,i, matrix_get(V,i,i)*b[0]);
    }

    A2  = matrix_complex_mult(A,A,1,NULL);
    A2n = matrix_complex_eye(dim, NULL);

    /* evaluate (10.33) */
    for(int i = 1; i <= M/2; i++)
    {
        /* A2n = A2n*A2 */
        matrix_complex_t *temp = matrix_complex_mult(A2n,A2,1,NULL);
        matrix_complex_free(A2n);
        A2n = temp;

        matrix_complex_add(U, A2n, b[2*i+1], NULL);
        matrix_complex_add(V, A2n, b[2*i],   NULL);
    }

    matrix_complex_free(A2);

    /* U = A*U */
    matrix_complex_mult(A,U,1,A2n);
    matrix_complex_free(U);
    U   = A2n;
    A2n = NULL;

    /* U = V+U, V = V-U */
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            const complex_t Uij = matrix_get(U,i,j);
            const complex_t Vij = matrix_get(V,i,j);

            matrix_set(U, i,j, Vij+Uij);
            matrix_set(V, i,j, Vij-Uij);
        }

    /* X = (V-U)^-1 * (V+U) */
    matrix_complex_invert(V);

    X = matrix_complex_mult(V,U,1,NULL);

    matrix_complex_free(U);
    matrix_complex_free(V);

    return X;
}


matrix_complex_t *_expm_ss(matrix_complex_t *A, const double norm)
{
    matrix_complex_t *A2 = NULL, *A4 = NULL, *A6 = NULL, *U = NULL, *V = NULL, *temp = NULL, *X = NULL;
    const int dim = A->rows;
    double const *b = b_list[13];
    const int s = ceil(log(norm/THETA_13)/M_LN2);

    matrix_complex_mult_scalar(A, pow(0.5,s));

    A2 = matrix_complex_mult(A, A, 1,NULL);
    A4 = matrix_complex_mult(A2,A2,1,NULL);
    A6 = matrix_complex_mult(A2,A4,1,NULL);

    /* calculate U */
    U = matrix_complex_alloc(dim,dim);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            matrix_set(U, i,j, b[13]*matrix_get(A6,i,j) + b[11]*matrix_get(A4,i,j) + b[9]*matrix_get(A2,i,j));

    temp = matrix_complex_mult(A6, U, 1, NULL);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            complex_t elem = matrix_get(temp,i,j);
            elem += b[7]*matrix_get(A6,i,j) + b[5]*matrix_get(A4,i,j) + b[3]*matrix_get(A2,i,j) + b[1]*(i == j ? 1 : 0);
            matrix_set(temp,i,j,elem);
        }
    matrix_complex_mult(A,temp,1,U);

    /* calculate V */
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            matrix_set(temp, i,j, b[12]*matrix_get(A6,i,j) + b[10]*matrix_get(A4,i,j) + b[8]*matrix_get(A2,i,j));

    V = matrix_complex_mult(A6, temp, 1, NULL);
    matrix_complex_free(temp);

    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            complex_t elem = matrix_get(V,i,j);
            elem += b[6]*matrix_get(A6,i,j) + b[4]*matrix_get(A4,i,j) + b[2]*matrix_get(A2,i,j) + b[0]*(i == j ? 1 : 0);
            matrix_set(V,i,j,elem);
        }

    matrix_complex_free(A2);
    matrix_complex_free(A4);
    matrix_complex_free(A6);


    /* U = V+U, V = V-U */
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            const complex_t Uij = matrix_get(U,i,j);
            const complex_t Vij = matrix_get(V,i,j);

            matrix_set(U, i,j, Vij+Uij);
            matrix_set(V, i,j, Vij-Uij);
        }

    /* X = (V-U)^-1 * (V+U) */
    matrix_complex_invert(V);

    X = matrix_complex_mult(V,U,1,NULL);

    matrix_complex_free(U);

    for(int i = 0; i < s; i++)
    {
        matrix_complex_t *temp;
        matrix_complex_mult(X,X,1,V);
        temp = X;
        X = V;
        V = temp;
    }

    matrix_complex_free(V);

    return X;
}

matrix_complex_t *matrix_complex_expm(matrix_complex_t *A)
{
    double norm;
    
    matrix_complex_norm(A, '1', &norm);

    if     (norm < THETA_3)
        return _expm_pade3579(A, 3);
    else if(norm < THETA_5)
        return _expm_pade3579(A, 5);
    else if(norm < THETA_7)
        return _expm_pade3579(A, 7);
    else if(norm < THETA_9)
        return _expm_pade3579(A, 9);
    else
        return _expm_ss(A,norm);
}

/** @}*/
