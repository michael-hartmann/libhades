/**
 * @file   expm.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   September, 2015
 * @brief  implement matrix exponential function
 */


#include <math.h>
#include <libhades.h>
#include <libhades/expm.h>

/* table 10.2 */
#define THETA_3  1.5e-2
#define THETA_5  2.5e-1
#define THETA_7  9.5e-1
#define THETA_9  2.1e0
#define THETA_13 5.4e0

#define return_error(cond,code) if((cond)) { ret = code; goto out; }

/* table 10.4 */
static double const b3[]  = { 120,60,12,1 };
static double const b5[]  = { 30240,15120,3360,420,30,1 };
static double const b7[]  = { 17297280, 8648640,1995840, 277200, 25200,1512, 56,1 };
static double const b9[]  = { 17643225600,8821612800,2075673600,302702400,30270240, 2162160,110880,3960,90,1 };
static double const b13[] = { 64764752532480000,32382376266240000,7771770303897600,1187353796428800,129060195264000,10559470521600,670442572800,33522128640,1323241920,40840800,960960,16380,182,1 };

/* lookup list for Pade coefficients */
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


/**
 * @brief Calculate matrix exponential of A using Pade approximation M
 *
 * This implements a part of algorithm 10.20.
 *
 * The calculation is performed inplace.
 *
 * @param [in,out] A square matrix
 * @param [in,out] M Pade approximation (must be 3,5,7 or 9)
 * @return 0 if successful
 * @return error otherwise
*/
static int _expm_pade3579(matrix_complex_t *A, int M)
{
    matrix_complex_t *U = NULL, *V = NULL, *A2 = NULL, *A2n = NULL, *workspace = NULL;
    int ret = 0;
    const int dim = A->rows;
    double const *b = b_list[M];

    /* evaluate (10.33) */
    {
        /* set U = b[1]*Id and V = b[0]*Id */
        U = matrix_complex_zeros(dim,dim, NULL);
        V = matrix_complex_zeros(dim,dim, NULL);
        return_error(U == NULL || V == NULL, LIBHADES_ERROR_OOM);
        for(int i = 0; i < dim; i++)
        {
            matrix_set(U, i,i, b[1]);
            matrix_set(V, i,i, b[0]);
        }

        /* A2 = A*A */
        A2 = matrix_complex_mult(A,A,1,NULL);
        return_error(A2 == NULL, LIBHADES_ERROR_OOM);

        A2n = matrix_complex_copy(A2,NULL);
        workspace = matrix_complex_alloc(dim,dim);
        return_error(A2 == NULL || workspace == NULL, LIBHADES_ERROR_OOM);

        matrix_complex_add(U, A2, b[3], NULL);
        matrix_complex_add(V, A2, b[2], NULL);
        for(int i = 2; i <= M/2; i++)
        {
            /* A2n = A2n*A2 */
            matrix_complex_mult(A2n,A2,1,workspace);
            matrix_complex_swap(A2n,workspace);

            matrix_complex_add(U, A2n, b[2*i+1], NULL);
            matrix_complex_add(V, A2n, b[2*i],   NULL);
        }

        matrix_complex_free(workspace);
        matrix_complex_free(A2);
        workspace = A2 = NULL;

        /* U = A*U */
        matrix_complex_mult(A,U,1,A2n);
        matrix_complex_free(U);
        U   = A2n;
        A2n = NULL;
    }

    /* U = V+U, V = V-U */
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            const complex_t Uij = matrix_get(U,i,j);
            const complex_t Vij = matrix_get(V,i,j);

            matrix_set(U, i,j, Vij+Uij);
            matrix_set(V, i,j, Vij-Uij);
        }

    /* A = (V-U)^-1 * (V+U) */
    ret = matrix_complex_invert(V); /* XXX solve linear system, not inv!!! */
    return_error(ret != 0, ret);
    matrix_complex_mult(V,U,1,A);

out:
    if(U != NULL)
        matrix_complex_free(U);
    if(V != NULL)
        matrix_complex_free(V);
    if(A2 != NULL)
        matrix_complex_free(A2);
    if(A2n != NULL)
        matrix_complex_free(A2n);
    if(workspace != NULL)
        matrix_complex_free(workspace);

    return ret;
}


/**
 * @brief Calculate A^(2^N)
 *
 * The calculation is performed inplace.
 *
 * @param [in,out] A complex square matrix
 * @param [in,out] work complex square matrix with same dimension as A
 * @param [in,out] N exponent
 */
static void _expm_square(matrix_complex_t *A, matrix_complex_t *work, const int N)
{
    for(int i = 0; i < N/2; i++)
    {
        matrix_complex_mult(A,   A,   1,work);
        matrix_complex_mult(work,work,1,A);
    }

    /* if N is odd, we have to do one more matrix multiplication */
    if((N % 2) == 1)
    {
        matrix_complex_mult(A,A,1,work);

        /* now the result is in work, so swap the matrices A and work */
        matrix_complex_swap(work,A);
    }
}


/** @brief Calculate matrix exponential using scaling and square method
 *
 * This function implements the scaling and square method using Pade
 * approximation.
 *
 * The calculation is performed inplace.
 *
 * @param [in,out] A square matrix
 * @param [in]     norm 1-norm of A
 * @return 0 if successful
 * @return error otherwise
 */
static int _expm_ss(matrix_complex_t *A, const double norm)
{
    matrix_complex_t *A2 = NULL, *A4 = NULL, *A6 = NULL, *U = NULL, *V = NULL, *temp = NULL;
    int ret = 0;
    const int dim = A->rows;
    double const *b = b_list[13];
    const int s = MAX(0,ceil(log(norm/THETA_13)/M_LN2));

    if(s != 0)
        matrix_complex_mult_scalar(A, pow(0.5,s));

    A2 = matrix_complex_mult(A, A, 1,NULL);
    return_error(A2 == NULL, LIBHADES_ERROR_OOM);

    A4 = matrix_complex_mult(A2,A2,1,NULL);
    return_error(A4 == NULL, LIBHADES_ERROR_OOM);

    A6 = matrix_complex_mult(A2,A4,1,NULL);
    return_error(A6 == NULL, LIBHADES_ERROR_OOM);

    /* calculate U */
    {
        U = matrix_complex_alloc(dim,dim);
        return_error(U == NULL, LIBHADES_ERROR_OOM);

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
                matrix_set(U, i,j, b[13]*matrix_get(A6,i,j) + b[11]*matrix_get(A4,i,j) + b[9]*matrix_get(A2,i,j));

        temp = matrix_complex_mult(A6, U, 1, NULL);
        return_error(temp == NULL, LIBHADES_ERROR_OOM);

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
            {
                complex_t elem = matrix_get(temp,i,j);
                elem += b[7]*matrix_get(A6,i,j) + b[5]*matrix_get(A4,i,j) + b[3]*matrix_get(A2,i,j) + b[1]*(i == j ? 1 : 0);
                matrix_set(temp,i,j,elem);
            }
        matrix_complex_mult(A,temp,1,U);
    }

    /* calculate V */
    {
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
                matrix_set(temp, i,j, b[12]*matrix_get(A6,i,j) + b[10]*matrix_get(A4,i,j) + b[8]*matrix_get(A2,i,j));

        V = matrix_complex_mult(A6, temp, 1, NULL);
        return_error(V == NULL, LIBHADES_ERROR_OOM);

        matrix_complex_free(temp);
        temp = NULL;

        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
            {
                complex_t elem = matrix_get(V,i,j);
                elem += b[6]*matrix_get(A6,i,j) + b[4]*matrix_get(A4,i,j) + b[2]*matrix_get(A2,i,j) + b[0]*(i == j ? 1 : 0);
                matrix_set(V,i,j,elem);
            }
    }

    matrix_complex_free(A2);
    matrix_complex_free(A4);
    matrix_complex_free(A6);
    A2 = A4 = A6 = NULL;

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
    ret = matrix_complex_invert(V); /* XXX solve linear system, not inv!!! */
    return_error(ret != 0, ret);

    matrix_complex_mult(V,U,1,A);
    matrix_complex_free(U);
    U = NULL;

    /* square */
    _expm_square(A, V, s);

out:
    if(A2 != NULL)
        matrix_complex_free(A2);
    if(A4 != NULL)
        matrix_complex_free(A4);
    if(A6 != NULL)
        matrix_complex_free(A6);
    if(U != NULL)
        matrix_complex_free(U);
    if(V != NULL)
        matrix_complex_free(V);
    if(temp != NULL)
        matrix_complex_free(temp);

    return ret;
}


/** @brief Calculate matrix exponential using scaling and square method and Pade approximation
 *
 * This function implements the scaling and square method using Pade
 * approximation.
 *
 * The calculation is performed inplace.
 *
 * @param [in,out] A square matrix
 * @return 0 if successful
 * @return error otherwise
 */
int matrix_complex_expm(matrix_complex_t *A)
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
