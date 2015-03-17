#ifndef LIBHADES_ODEINT_H
#define LIBHADES_ODEINT_H

typedef struct {
    void (*f)(matrix_t *ft, matrix_t *y, double t, void *args);
    void (*Jf)(matrix_t *, matrix_t *, double, void *);
    void *args;

    int Jf_of_t;

    double epsilon;
    double maxiter;
} rk4_symplectic_t;


int rk4        (void (*f)(double t, matrix_t         *y, matrix_t         *ft, void *args), matrix_t         *yn, double t0, double t, double h, void *args);
int rk4_complex(void (*f)(double t, matrix_complex_t *y, matrix_complex_t *ft, void *args), matrix_complex_t *yn, double t0, double t, double h, void *args);

void rk4_symplectic_init(rk4_symplectic_t *self, void (*f)(matrix_t *, matrix_t *, double, void *), void (*Jf)(matrix_t *, matrix_t *, double, void *), void *args);
void rk4_symplectic_free(rk4_symplectic_t *self);

void rk4_symplectic_set_maxiter(rk4_symplectic_t *self, int maxiter);
void rk4_symplectic_set_epsilon(rk4_symplectic_t *self, double epsilon);
void rk4_symplectic_set_Jf_of_t(rk4_symplectic_t *self, int boolean);

int rk4_symplectic_integrate(rk4_symplectic_t *self, matrix_t *yn, double t, double t0, int steps);

#endif
