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

/** \defgroup odeint Integration of ordinary differential equations
 *  @{
 */

#define RK4_INIT(FUNCTION_NAME, RK4_TYPE, RK4_FREE, MATRIX_TYPE, MATRIX_ALLOC) \
    int FUNCTION_NAME(RK4_TYPE *self, void (*f)(MATRIX_TYPE *, MATRIX_TYPE *, double, void *), int dim, void *args) \
    { \
        self->f    = f; \
        self->args = args; \
\
        for(int i = 0; i < sizeof(self->workspace)/sizeof(self->workspace[0]); i++) \
            self->workspace[i] = NULL; \
\
        for(int i = 0; i < sizeof(self->workspace)/sizeof(self->workspace[0]); i++) \
        { \
            self->workspace[i] = MATRIX_ALLOC(dim,1); \
            if(self->workspace[i] == NULL) \
            { \
                RK4_FREE(self); \
                return LIBHADES_ERROR_OOM; \
            } \
        } \
\
        return 0; \
    }


/** @brief Initialize workspace for explicit 4th order Runge Kutta method
 *
 * This function will initialize the Runge Kutta 4th order integration to solve
 * the differential equation
 *      d/dt y = f(y,t),
 * where y is a vector and f is a arbitrary function for initial condition
 * y(t0) = y0. The differential equation is solved using explicit Runge Kutta
 * 4th order.
 *
 * The callback function f must have type
 *      (*f)(matrix_complex_t *dest, matrix_complex_t *src, double t, void *args).
 * The result must be stored in dest, the input vector is src, the time is t,
 * and args is an arbitrary pointer given to the function.
 *
 * @param [in,out] self explicit Runge-Kutta 4th order object for real differential equation
 * @param [in] f callback for rhs of differential equation
 * @param [in] dim dimension of differential equation (rows of vector y)
 * @param [in] args arbitrary pointer that will be passed to f
 *
 * @return ret 0 if successfull, otherwise LIBHADES_ERROR_OOM
 */
RK4_INIT(rk4_init, rk4_t, rk4_free, matrix_t, matrix_alloc)


/** @brief Initialize workspace for explicit 4th order Runge Kutta method
 *
 * See rk4_init
 */
RK4_INIT(rk4_complex_init, rk4_complex_t, rk4_complex_free, matrix_complex_t, matrix_complex_alloc)

#define RK4_FREE(FUNCTION_NAME, RK4_TYPE, MATRIX_FREE) \
void FUNCTION_NAME(RK4_TYPE *self) \
{ \
    for(int i = 0; i < 5; i++) \
        if(self->workspace[i] != NULL) \
            MATRIX_FREE(self->workspace[i]); \
}

/** @brief Free memory allocated by rk4_init
 *
 * @param self rk4_t object
 **/
RK4_FREE(rk4_free, rk4_t, matrix_free)

/** @brief Free memory allocated by rk4_complex_init
 *
 * @param self rk4_t object
 **/
RK4_FREE(rk4_complex_free, rk4_complex_t, matrix_complex_free)


#define RK4_INTEGRATE(FUNCTION_NAME, RK4_TYPE, TYPE, MATRIX_TYPE, MATRIX_ALLOC, MATRIX_ADD, MATRIX_FREE) \
int FUNCTION_NAME(RK4_TYPE *self, MATRIX_TYPE *y0, double t, double t0, int steps) \
{ \
    const int N = y0->rows; \
    double tn = t0; \
    const double h = (t-t0)/steps;\
    MATRIX_TYPE *yn = y0; \
\
    MATRIX_TYPE *k1 = self->workspace[0]; \
    MATRIX_TYPE *k2 = self->workspace[1]; \
    MATRIX_TYPE *k3 = self->workspace[2]; \
    MATRIX_TYPE *k4 = self->workspace[3]; \
    MATRIX_TYPE *z  = self->workspace[4]; \
\
    TYPE *ynM = yn->M; \
    TYPE *k1M = k1->M; \
    TYPE *k2M = k2->M; \
    TYPE *k3M = k3->M; \
    TYPE *k4M = k4->M; \
\
    for(int n = 0; n < steps; n++) \
    { \
        tn = t0+n*h; \
\
        /* k1 */ \
        self->f(k1, yn, tn+h, self->args); \
\
        /* k2 */ \
        MATRIX_ADD(yn, k1, h/2, z); \
        self->f(k2, z, tn+h/2, self->args); \
\
        /* k3 */ \
        MATRIX_ADD(yn, k2, h/2, z); \
        self->f(k3, z, tn+h/2, self->args); \
\
        /* k4 */ \
        MATRIX_ADD(yn, k3, h, z); \
        self->f(k4, z, tn+h, self->args); \
\
        for(int j = 0; j < N; j++) \
            ynM[j] += h/6*(k1M[j]+2*(k2M[j]+k3M[j])+k4M[j]); \
    } \
\
    return 0; \
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
RK4_INTEGRATE(rk4_integrate, rk4_t, double, matrix_t, matrix_alloc, matrix_add, matrix_free);


/** @brief Ruge Kutta 4th order for complex ODE
 *
 * See rk4.
 */
RK4_INTEGRATE(rk4_complex_integrate, rk4_complex_t, complex_t, matrix_complex_t, matrix_complex_alloc, matrix_complex_add, matrix_complex_free);


#define RK4_INTEGRATE_VARIABLE_STEPSIZE(FUNCTION_NAME, RK4_TYPE, RK4_INTEGRATE_FUNC, MATRIX_TYPE) \
    int FUNCTION_NAME(RK4_TYPE *self, MATRIX_TYPE *yn, double t, double t0, double (*adapt)(void (*f)(MATRIX_TYPE *, MATRIX_TYPE *, double, void *), MATRIX_TYPE *y, double t, void *args)) \
    { \
        int done = 0, steps = 0; \
        double tn = t0; \
\
        while(!done) \
        { \
            double h = adapt(self->f, yn, tn, self->args); \
            if((tn+h) > t) \
            { \
                h = t-tn; \
                done = 1; \
            } \
\
            RK4_INTEGRATE_FUNC(self, yn, tn+h, tn, 1); \
            tn += h; \
            steps++; \
        } \
\
        return steps; \
    }


RK4_INTEGRATE_VARIABLE_STEPSIZE(rk4_integrate_variable_stepsize, rk4_t, rk4_integrate, matrix_t)
RK4_INTEGRATE_VARIABLE_STEPSIZE(rk4_complex_integrate_variable_stepsize, rk4_complex_t, rk4_complex_integrate, matrix_complex_t)


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
int rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args)
{
    return rk4_symplectic_init_full(self, f, dim, args, 50000, 1e-16, NULL);
}

int rk4_symplectic_init_full(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args, int maxiter, double epsilon, double (adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args))
{
    self->f    = f;
    self->args = args;

    self->maxiter = maxiter;
    self->epsilon = epsilon;

    self->adapt = adapt;

    /* workspace */
    for(int i = 0; i < 8; i++)
        self->workspace[i] = NULL;

    for(int i = 0; i < 8; i++)
    {
        self->workspace[i] = matrix_alloc(dim,1);
        if(self->workspace[i] == NULL)
        {
            rk4_symplectic_free(self);
            return LIBHADES_ERROR_OOM;
        }
    }

    return 0;
}


/** @brief Free memory of object */
void rk4_symplectic_free(rk4_symplectic_t *self)
{
    for(int i = 0; i < 8; i++)
        if(self->workspace[i] != NULL)
            matrix_free(self->workspace[i]);
}


static inline void interpolate(const double t1, const double t2, const double x[3], const double f[3], double *w1, double *w2)
{
    double a0 = f[0];
    double a1 = (f[1]-a0)/( x[1]-x[0] );
    double a2 = ( (f[2]-f[0])/(x[2]-x[0]) - a1 )/(x[2]-x[1]);

    *w1 = a0 + (t1-x[0])*( a1 + a2*(t1-x[1]) );
    *w2 = a0 + (t2-x[0])*( a1 + a2*(t2-x[1]) );
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
    const int dim = yn->rows;
    const double epsilon = self->epsilon;
    const double maxiter = self->maxiter;
    double tn = t0;
    double h = (t-t0)/steps;
    int bye = 0;

    matrix_t *Z1      = self->workspace[0];
    matrix_t *Z1_last = self->workspace[1];
    matrix_t *Z2      = self->workspace[2];
    matrix_t *Z2_last = self->workspace[3];

    matrix_t *ynpZ1   = self->workspace[4];
    matrix_t *ynpZ2   = self->workspace[5];

    matrix_t *Y1      = self->workspace[6];
    matrix_t *Y2      = self->workspace[7];

    /* init */
    matrix_setall(Z1, 0);
    matrix_setall(Z2, 0);

    if(self->adapt != NULL)
        h = self->adapt(self->f, yn, t, self->args);

    if((t-tn) < h)
    {
        h = t-tn;
        bye = 1;
    }

    while(1)
    {
        /* fixed point iteration */
        for(int j = 0; ; j++)
        {
            double norm_Z1, norm_Z2;

            Z1_last = self->workspace[j     % 2];
            Z1      = self->workspace[(j+1) % 2];
            Z2_last = self->workspace[j     % 2 + 2];
            Z2      = self->workspace[(j+1) % 2 + 2];

            /* ynpZ1 = yn+Z1, ynpZ2 = yn+Z2 */
            matrix_add(yn, Z1_last, 1, ynpZ1);
            matrix_add(yn, Z2_last, 1, ynpZ2);

            /* Z1 = a11*h*f(tn+c1*h, yn+Z1) + a12*h*f(tn+c1*h, yn+Z2) */
            self->f(Z1, ynpZ1, tn+c1*h, self->args);
            matrix_mult_scalar(Z1, a11*h);
            self->f(Y1, ynpZ2, tn+c1*h, self->args);
            matrix_add(Z1, Y1, a12*h, NULL);

            /* Z2 = a21*h*f(tn+c2*h, yn+Z1) + a22*h*f(tn+c2*h, yn+Z2) */
            self->f(Z2, ynpZ1, tn+c2*h, self->args);
            matrix_mult_scalar(Z2, a21*h);
            self->f(Y1, ynpZ2, tn+c2*h, self->args);
            matrix_add(Z2, Y1, a22*h, NULL);

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
                fprintf(stderr, "Convergence error, h=%g\n", h);
                return j;
            }
        }

        matrix_add(yn, Z1, 1, ynpZ1);
        matrix_add(yn, Z2, 1, ynpZ2);

        self->f(Y1, ynpZ1, tn+c1*h, self->args);
        self->f(Y2, ynpZ2, tn+c2*h, self->args);

        matrix_add(yn, Y1, h*b1, NULL);
        matrix_add(yn, Y2, h*b2, NULL);

        if(bye)
            return 0;

        tn += h;

        if(self->adapt != NULL)
            h = self->adapt(self->f, yn, t, self->args);

        if((t-tn) < h)
        {
            h = t-tn;
            bye = 1;
        }

        /* determine starting approximations */
        {
            const double t[] = { tn, tn+h*c1, tn+h*c2 };
            double f[3];
            f[0] = 0;

            for(int k = 0; k < dim; k++)
            {
                double w1,w2;
                f[1] = matrix_get(Z1, k,0);
                f[2] = matrix_get(Z2, k,0);

                interpolate(tn+c1*h, tn+c2*h, t, f, &w1, &w2);

                matrix_set(Z1, k,0, w1);
                matrix_set(Z2, k,0, w2);
            }
        }
    }
}

/** @}*/
