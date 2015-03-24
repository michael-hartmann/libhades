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

#define b 0.5 /* b1 = b1 = b */

#define c1 0.21132486540518713
#define c2 0.7886751345948129

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
 * @param self [out] symplectic Runge Kutta object
 * @param f    [in]  callback of function f
 * @param args [in]  arbitrary data that is passed to f and Jf
 */
void rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), void *args)
{
    rk4_symplectic_init_full(self, f, args, 50000, 1e-12, NULL);
}

void rk4_symplectic_init_full(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), void *args, int maxiter, double epsilon, double (adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args))
{
    self->f    = f;
    self->args = args;

    self->maxiter = maxiter;
    self->epsilon = 1e-12;

    self->adapt = adapt;
}


/** @brief Free memory of object */
void rk4_symplectic_free(rk4_symplectic_t *self)
{
    return;
}


/** @brief Integrate differential equation
 *
 * 
 * @param self  [in,out] symplectic Runge Kutta object
 * @param yn    [in,out] at start initial conditions yn=y(t0), at the end yn=y(t)
 * @param t     [in]     end time
 * @param t0    [in]     initial time
 * @param steps [in]     number of steps used
 * @return success 0 if successful, LIBHADES_ERROR_OOM if out of memory, >0 if convergence failed
 */
int rk4_symplectic_integrate(rk4_symplectic_t *self, matrix_t *yn, double t, double t0, int steps)
{
    int ret = 0;
    const int rows = yn->rows;
    const double epsilon = self->epsilon;
    const double maxiter = self->maxiter;
    double tn = t0;

    matrix_t *Z1      = matrix_alloc(rows,1);
    matrix_t *Z2      = matrix_alloc(rows,1);
    matrix_t *Z1_last = matrix_alloc(rows,1);
    matrix_t *Z2_last = matrix_alloc(rows,1);

    matrix_t *ynpZ1 = matrix_alloc(rows, 1);
    matrix_t *ynpZ2 = matrix_alloc(rows, 1);

    matrix_t *f1 = matrix_alloc(rows, 1);
    matrix_t *f2 = matrix_alloc(rows, 1);

    if(Z1 == NULL || Z2 == NULL || ynpZ1 == NULL || ynpZ2 == NULL || f1 == NULL || f2 == NULL)
    {
        ret = LIBHADES_ERROR_OOM;
        goto out;
    }

    while(1)
    {
        double h;
        int bye = 0;

        if(self->adapt == NULL)
            h = (t-t0)/steps;
        else
            h = self->adapt(self->f, yn, t, self->args);

        if((t-tn) < h)
        {
            h = t-tn;
            bye = 1;
        }

        /* fixed point iteration */
        /* XXX find better starting values */
        matrix_setall(Z1, 0);
        matrix_setall(Z2, 0);

        for(int j = 0; ; j++)
        {
            double norm_Z1, norm_Z2;

            matrix_copy(Z1, Z1_last);
            matrix_copy(Z2, Z2_last);

            matrix_add(yn, Z1, 1, ynpZ1);
            matrix_add(yn, Z2, 1, ynpZ2);

            self->f(Z1, ynpZ1, tn+c1*h, self->args);
            matrix_mult_scalar(Z1, a11*h);
            self->f(f1, ynpZ2, tn+c1*h, self->args);
            matrix_add(Z1, f1, a12*h, NULL);

            self->f(Z2, ynpZ1, tn+c2*h, self->args);
            matrix_mult_scalar(Z2, a21*h);
            self->f(f1, ynpZ2, tn+c2*h, self->args);
            matrix_add(Z2, f1, a22*h, NULL);

            matrix_add(Z1_last, Z1, -1, NULL);
            matrix_norm(Z1_last, 'F', &norm_Z1);

            if(norm_Z1 < epsilon)
            {
                matrix_add(Z2_last, Z2, -1, NULL);
                matrix_norm(Z2_last, 'F', &norm_Z2);

                if(norm_Z2 < epsilon)
                    break;
            }
            if(j >= maxiter)
            {
                fprintf(stderr, "Convergence error\n");
                ret = j;
                goto out;
            }
        }

        matrix_add(yn, Z1, 1, ynpZ1);
        matrix_add(yn, Z2, 1, ynpZ2);

        self->f(f1, ynpZ1, tn+c1*h, self->args);
        self->f(f2, ynpZ2, tn+c2*h, self->args);

        matrix_add(f1, f2, 1, NULL);
        matrix_add(yn, f1, h*b, NULL);

        tn += h;

        if(bye)
            break;
    }

out:
    if(Z1 != NULL)
        matrix_free(Z1);
    if(Z2 != NULL)
        matrix_free(Z2);

    if(ynpZ1 != NULL)
        matrix_free(ynpZ1);
    if(ynpZ2 != NULL)
        matrix_free(ynpZ2);

    if(f1 != NULL)
        matrix_free(f1);
    if(f2 != NULL)
        matrix_free(f2);

    return ret;
}

/** @}*/
