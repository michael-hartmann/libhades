#ifndef LIBHADES_OPTIMIZE_H
#define LIBHADES_OPTIMIZE_H

int newton_mdim(void (*f)(matrix_t *, matrix_t *, void *), void (*Jacobian)(matrix_t *, matrix_t *, void *), matrix_t *xn, double eps, int maxiter, void *args);

#endif
