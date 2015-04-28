#ifndef LIBHADES_ODEINT_H
#define LIBHADES_ODEINT_H

#include <libhades.h>

typedef struct {
    void (*f)(matrix_t *ft, matrix_t *y, double t, void *args);
    int maxiter;
    double epsilon;
    void *args;
    double (*adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args);

    /* workspace */
    matrix_t *workspace[8];
} rk4_symplectic_t;

typedef struct {
    void (*f)(matrix_t *ft, matrix_t *y, double t, void *args);
    int maxiter;
    double epsilon;
    void *args;
    double (*adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args);

    /* workspace */
    matrix_t *workspace[5];
} rk4_t;

typedef struct {
    void (*f)(matrix_complex_t *ft, matrix_complex_t *y, double t, void *args);
    int maxiter;
    double epsilon;
    void *args;
    double (*adapt)(void (*f)(matrix_complex_t *, matrix_complex_t *, double, void *), matrix_complex_t *y, double t, void *args);

    /* workspace */
    matrix_complex_t *workspace[5];
} rk4_complex_t;

int rk4_init        (rk4_t         *self, void (*f)(matrix_t *, matrix_t *, double, void *),                 int dim, void *args);
int rk4_complex_init(rk4_complex_t *self, void (*f)(matrix_complex_t *, matrix_complex_t *, double, void *), int dim, void *args);

void rk4_free        (rk4_t         *self);
void rk4_complex_free(rk4_complex_t *self);

int rk4_integrate        (rk4_t         *self, matrix_t         *yn, double t, double t0, int steps);
int rk4_complex_integrate(rk4_complex_t *self, matrix_complex_t *yn, double t, double t0, int steps);

int rk4_integrate_variable_stepsize        (rk4_t         *self, matrix_t         *yn, double t, double t0, double (*adapt)(void (*f)(matrix_t         *, matrix_t         *, double, void *), matrix_t         *y, double t, void *args));
int rk4_complex_integrate_variable_stepsize(rk4_complex_t *self, matrix_complex_t *yn, double t, double t0, double (*adapt)(void (*f)(matrix_complex_t *, matrix_complex_t *, double, void *), matrix_complex_t *y, double t, void *args));


int rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args);
int rk4_symplectic_init_full(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args, int maxiter, double epsilon, double (adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args));

void rk4_symplectic_free(rk4_symplectic_t *self);

int rk4_symplectic_integrate(rk4_symplectic_t *self, matrix_t *yn, double t, double t0, int steps);

#endif
