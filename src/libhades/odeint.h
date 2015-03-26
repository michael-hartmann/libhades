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


int rk4        (void (*f)(double t, matrix_t         *y, matrix_t         *ft, void *args), matrix_t         *yn, double t0, double t, double h, void *args);
int rk4_complex(void (*f)(double t, matrix_complex_t *y, matrix_complex_t *ft, void *args), matrix_complex_t *yn, double t0, double t, double h, void *args);

int rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args);
int rk4_symplectic_init_full(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), int dim, void *args, int maxiter, double epsilon, double (adapt)(void (*f)(matrix_t *, matrix_t *, double, void *), matrix_t *y, double t, void *args));

void rk4_symplectic_free(rk4_symplectic_t *self);

int rk4_symplectic_integrate(rk4_symplectic_t *self, matrix_t *yn, double t, double t0, int steps);

#endif
