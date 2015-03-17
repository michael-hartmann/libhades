/**
 * @file   odeint.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   March, 2015
 * @brief  ordinary and symplectic Runge Kutta 4th order
 */


#include <math.h>
#include <libhades.h>
#include <libhades/optimize.h>
#include <libhades/odeint.h>

/* Butcher tableau */
#define a11 0.25
#define a12 -0.038675134594812866
#define a21 +0.5386751345948129
#define a22 0.25

#define b1 0.5
#define b2 0.5

#define c1 0.21132486540518713
#define c2 0.7886751345948129


typedef struct {
    double tn,h;
    void (*f)(matrix_t *ft, matrix_t *y, double t, void *args);
    void (*Jf)(matrix_t *, matrix_t *, double, void *);
    matrix_t *yn;
    void *args;
    int Jf_of_t;

    matrix_t *f_xi1, *f_xi2, *f11, *f12, *f21, *f22;
} rk4_params;

/** \defgroup odeint Integration of ordinary differential equations
 *  @{
 */

/* d/dt y = f(t,y), y(0)=y0 */
#define RK4(FUNCTION_NAME, TYPE, MATRIX_TYPE, MATRIX_ALLOC, MATRIX_ADD, MATRIX_FREE) \
int FUNCTION_NAME(void (*f)(double t, MATRIX_TYPE *y, MATRIX_TYPE *ft, void *args), MATRIX_TYPE *y0, double t0, double t, double h, void *args) \
{ \
    const int N = y0->rows; \
    int ret = 0; \
    double tn = t0; \
    const int steps = ceil((t-t0)/h); \
    h = t/steps; \
    MATRIX_TYPE *yn = y0; \
\
    MATRIX_TYPE *k1 = MATRIX_ALLOC(N,1); \
    MATRIX_TYPE *k2 = MATRIX_ALLOC(N,1); \
    MATRIX_TYPE *k3 = MATRIX_ALLOC(N,1); \
    MATRIX_TYPE *k4 = MATRIX_ALLOC(N,1); \
\
    TYPE *ynM = yn->M; \
    TYPE *k1M = k1->M; \
    TYPE *k2M = k2->M; \
    TYPE *k3M = k3->M; \
    TYPE *k4M = k4->M; \
\
    MATRIX_TYPE *z = MATRIX_ALLOC(N,1); \
\
    if(k1 == NULL || k2 == NULL || k3 == NULL || k4 == NULL || z == NULL)\
    { \
        ret = LIBHADES_ERROR_OOM; \
        goto out; \
    } \
\
    for(int n = 0; n < steps; n++) \
    { \
        tn = t0+n*h; \
\
        /* k1 */ \
        f(tn+h,   yn, k1, args); \
\
        /* k2 */ \
        MATRIX_ADD(yn, k1, h/2, z); \
        f(tn+h/2, z, k2, args); \
\
        /* k3 */ \
        MATRIX_ADD(yn, k2, h/2, z); \
        f(tn+h/2, z, k3, args); \
\
        /* k4 */ \
        MATRIX_ADD(yn, k3, h, z); \
        f(tn+h,   z, k4, args); \
\
        for(int j = 0; j < N; j++) \
            ynM[j] += h/6*(k1M[j]+2*(k2M[j]+k3M[j])+k4M[j]); \
    } \
 \
    out: \
    MATRIX_FREE(k1); \
    MATRIX_FREE(k2); \
    MATRIX_FREE(k3); \
    MATRIX_FREE(k4); \
    MATRIX_FREE(z); \
\
    return ret; \
} \


/** @brief Ruge Kutta 4th order for real ODE
 *
 * Integrate vector yn=y(t0) from t=t0 to t for differential equation:
 * dy/dt = f(y,t)
 *
 * f is a callback function with arguments
 *   - t:    time
 *   - y:    vector y
 *   - ft:   this vector will hold f(y,t)
 *   - args: void pointer given to function
 *
 * @param [in] f callback function
 * @param [in,out] yn y(t0) at start, y(t) afterwards
 * @param [in] t0 initial time
 * @param [in] t end time
 * @param [in] h precision
 * @param [in] args pointer given to function f
 */
RK4(rk4, double, matrix_t, matrix_alloc, matrix_add, matrix_free);

/** @brief Ruge Kutta 4th order for complex ODE
 *
 * See rk4.
 */
RK4(rk4_complex, complex_t, matrix_complex_t, matrix_complex_alloc, matrix_complex_add, matrix_complex_free);


/** function that calculates Jacobian of the non-linear system of equations */
static void rk4_symplectic_J(matrix_t *dest, matrix_t *xi, void *p)
{
    const int dim = xi->rows;
    const int dimby2 = dim/2;
    matrix_t *f11, *f12, *f21, *f22;

    rk4_params *params = p;
    const double tn = params->tn;
    const double h  = params->h;
    f11 = params->f11;
    f12 = params->f12;
    f21 = params->f21;
    f22 = params->f22;

    matrix_t xi1 = {
        .rows    = xi->rows/2,
        .columns = 1,
        .min     = xi->rows/2, 
        .size    = xi->rows/2,
        .type    = 0,
        .view    = 1,
        .M       = xi->M
    };

    matrix_t xi2 = {
        .rows    = xi->rows/2,
        .columns = 1,
        .min     = xi->rows/2, 
        .size    = xi->rows/2,
        .type    = 0,
        .view    = 1,
        .M       = xi->M+xi->rows/2
    };


    params->Jf(f11, &xi1, tn+c1*h, params->args);
    params->Jf(f12, &xi2, tn+c1*h, params->args);

    if(params->Jf_of_t)
    {
        params->Jf(f21, &xi1, tn+c2*h, params->args);
        params->Jf(f22, &xi2, tn+c2*h, params->args);
    }
    else
    {
        f21 = f11;
        f22 = f12;
    }

    for(int i = 0; i < dimby2; i++)
    {
        for(int j = 0; j < dimby2; j++)
        {
            matrix_set(dest, i,j,               h*a11*matrix_get(f11,i,j));
            matrix_set(dest, i,j+dimby2,        h*a12*matrix_get(f12,i,j));
            matrix_set(dest, i+dimby2,j,        h*a21*matrix_get(f21,i,j));
            matrix_set(dest, i+dimby2,j+dimby2, h*a22*matrix_get(f22,i,j));
        }

        matrix_set(dest, i,i, matrix_get(dest,i,i)-1);
        matrix_set(dest, i+dimby2,i+dimby2, matrix_get(dest,i+dimby2,i+dimby2)-1);
    }
}


/** function that calculates function of the non-linear system of equations */
static void rk4_symplectic_f(matrix_t *dest, matrix_t *xi, void *p)
{
    const int rows = xi->rows;
    int i;
    matrix_t *f_xi1, *f_xi2;

    rk4_params *params = p;
    double tn = params->tn;
    double h = params->h;
    matrix_t *yn = params->yn;

    f_xi1 = params->f_xi1;
    f_xi2 = params->f_xi2;

    matrix_t xi_j = {
        .rows    = rows/2,
        .columns = 1,
        .min     = rows/2, 
        .size    = rows/2,
        .type    = 0,
        .view    = 1,
        .M       = NULL
    };

    for(i = 0; i < rows; i++)
        matrix_set(dest, i,0, matrix_get(yn, i % (rows/2),0)-matrix_get(xi, i,0));
    
    xi_j.M = xi->M;
    params->f(f_xi1, &xi_j, tn+c1*h, params->args);
    xi_j.M = xi->M+rows/2;
    params->f(f_xi2, &xi_j, tn+c1*h, params->args);

    for(i = 0; i < rows/2; i++)
        matrix_set(dest, i,0, matrix_get(dest,i,0)+h*( a11*matrix_get(f_xi1, i,0) + a12*matrix_get(f_xi2, i,0) ));

    xi_j.M = xi->M;
    params->f(f_xi1, &xi_j, tn+c2*h, params->args);
    xi_j.M = xi->M+rows/2;
    params->f(f_xi2, &xi_j, tn+c2*h, params->args);

    for(i = 0; i < rows/2; i++)
        matrix_set(dest, i+rows/2,0, matrix_get(dest,i+rows/2,0)+h*( a21*matrix_get(f_xi1, i,0) + a22*matrix_get(f_xi2, i,0) ));
}


/** @brief Initialize symplectic Runge Kutta 4th order integration
 *
 * This function will initialize the Runge Kutta 4th order integration to solve
 * the the differential equation d/dt y = f(y,t), where y is a vector and f is
 * a arbitrary function. The differential equation is solved using implicit,
 * symplectic Runge Kutta 4th order.
 *
 * The function f must be of type
 *      void f(matrix_t *ft, matrix_t *y, double t, void *args).
 * The result of f will be stored in ft, y and t are the arguments, args is a
 * pointer to arbitrary data.
 *
 * In order to integrate the differential equation, the Jacobian of the
 * function f is needed. The function Jf must be of type
 *      void Jf(matrix_t *Jft, matrix_t *y, double t, void *args).
 * The result of Jf(y,t) is stored in Jft, y and t are the arguments, and args
 * is a pointer to arbitrary data.
 *
 * @param self [out] symplectic Runge Kutta object
 * @param f    [in]  callback of function f
 * @param Jf   [in]  callback of Jacobian of f
 * @param args [in]  arbitrary data that is passed to f and Jf
 */
void rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), void (*Jf)(matrix_t *, matrix_t *, double, void *), void *args)
{
    self->f    = f;
    self->Jf   = Jf;
    self->args = args;

    self->Jf_of_t = 1;

    self->maxiter = 50;
    self->epsilon = 1e-12;
}


/** @brief Free memory of object */
void rk4_symplectic_free(rk4_symplectic_t *self)
{
    return;
}


/** @brief Set maximum number of iterations for Newton's method
 *
 * This is the maximum number of iterations for the Newton method to solve the
 * (in general) non linear system of equations. Default: 50
 *
 * @param self [in,out] symplectic Runge Kutta object
 * @param maxiter [in]  maximum number of iteration
 */
void rk4_symplectic_set_maxiter(rk4_symplectic_t *self, int maxiter)
{
    self->maxiter = maxiter;
}


/** @brief Set epsilon for Newton's method
 *
 * The iteration of Newton's method will stop if ||x_n+1 - x_n|| < epsilon.
 * The norm used is the Frobenius norm.
 *
 * By default: epsilon = 1e-12
 *
 * @param self [in,out] symplectic Runge Kutta object
 * @param epsilon [in]  epsilon
 */
void rk4_symplectic_set_epsilon(rk4_symplectic_t *self, double epsilon)
{
    self->epsilon = epsilon;
}


/** @brief Assume Jf is / is not time dependent
 *
 * If Jt is not a function of t, you may set:
 *      rk4_symplectic_set_Jf_of_t(self, 0)
 * This will increase performance as the Jacobian needs to be evaluated less
 * often.
 *
 * By default: assume Jt depends on time
 *
 * @param self [in,out] symplectic Runge Kutta object
 * @param boolean [in]  boolean
 */
void rk4_symplectic_set_Jf_of_t(rk4_symplectic_t *self, int boolean)
{
    self->Jf_of_t = boolean;
}


/** @brief Integrate differential equation
 *
 * 
 * @param self  [in,out] symplectic Runge Kutta object
 * @param yn    [in,out] at start initial conditions yn=y(t0), at the end yn=y(t)
 * @param t     [in]     end time
 * @param t0    [in]     initial time
 * @param steps [in]     number of steps used
 * @return success 0 if successful, LIBHADES_ERROR_OOM if out of memory, >0 if Newton's method failed
 */
int rk4_symplectic_integrate(rk4_symplectic_t *self, matrix_t *yn, double t, double t0, int steps)
{
    int ret = 0;
    const int rows = yn->rows;
    const double h = (t-t0)/steps;
    matrix_t *xi    = matrix_alloc(2*rows,1);
    matrix_t *f1    = matrix_alloc(rows, 1);
    matrix_t *f2    = matrix_alloc(rows, 1);

    matrix_t *f_xi1 = matrix_alloc(rows, 1);
    matrix_t *f_xi2 = matrix_alloc(rows, 1);

    matrix_t *f11 = matrix_alloc(2*rows, 2*rows);
    matrix_t *f12 = matrix_alloc(2*rows, 2*rows);
    matrix_t *f21 = matrix_alloc(2*rows, 2*rows);
    matrix_t *f22 = matrix_alloc(2*rows, 2*rows);

    if(xi == NULL || f1 == NULL || f2 == NULL || f_xi1 == NULL || f_xi2 == NULL || f11 == NULL || f12 == NULL || f21 == NULL || f22 == NULL)
    {
        ret = LIBHADES_ERROR_OOM;
        goto out;
    }

    matrix_t xi_j = {
        .rows    = rows,
        .columns = 1,
        .min     = 1,
        .size    = rows,
        .type    = 0,
        .view    = 1,
        .M       = NULL,
    };


    rk4_params p = {
        .tn   = 0,
        .h    = h,
        .f    = self->f,
        .yn   = yn,
        .Jf   = self->Jf,
        .args = self->args,

        .Jf_of_t = self->Jf_of_t,

        .f_xi1   = f_xi1,
        .f_xi2   = f_xi2,

        .f11 = f11,
        .f12 = f12,
        .f21 = f21,
        .f22 = f22
    };


    for(int i = 0; i < steps; i++)
    {
        const double tn = t0+i*h;
        for(int j = 0; j < rows; j++)
        {
            const double elem = matrix_get(yn,j,0);
            matrix_set(xi, j,     0, elem);
            matrix_set(xi, j+rows,0, elem);
        }

        p.tn = tn;

        ret = newton_mdim(rk4_symplectic_f, rk4_symplectic_J, xi, self->epsilon, self->maxiter, &p);
        if(ret < 0)
            goto out;

        xi_j.M = xi->M;
        self->f(f1, &xi_j, tn+c1*h, self->args);
        matrix_add(yn, f1, b1*h, NULL);

        xi_j.M = xi->M+rows;
        self->f(f2, &xi_j, tn+c2*h, self->args);
        matrix_add(yn, f2, b2*h, NULL);
    }

out:
    if(xi != NULL)
        matrix_free(xi);
    if(f1 != NULL)
        matrix_free(f1);
    if(f2 != NULL)
        matrix_free(f2);

    if(f_xi1 != NULL)
        matrix_free(f_xi1);
    if(f_xi2 != NULL)
        matrix_free(f_xi2);

    if(f11 != NULL)
        matrix_free(f11);
    if(f12 != NULL)
        matrix_free(f12);
    if(f21 != NULL)
        matrix_free(f21);
    if(f22 != NULL)
        matrix_free(f22);

    return ret;
}

/** @}*/
