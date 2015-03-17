#include <libhades.h>

/** @brief Find roots of a multi-dimensional function f(x)
 *
 * @param [in] f                function f(x)
 * @param [in] Jacobian         Jacobian of function f(x)
 * @param [in,out] xn           on input the start for the iteration, on output the estimate for the root
 * @param [in,out] eps          eps
 * @param [in,out] maxiter      maximum number of iterations
 * @param [in,out] args         pointer that will be given to f and Jacobian
 *
 * @retval 0 if successful
 */
int newton_mdim(void (*f)(matrix_t *, matrix_t *, void *), void (*Jacobian)(matrix_t *, matrix_t *, void *), matrix_t *xn, double eps, int maxiter, void *args)
{
    const int dim = xn->rows;
    int i, converged = 0;

    matrix_t *J = matrix_alloc(dim,dim);
    matrix_t *b = matrix_alloc(dim,1);

    for(i = 0; i < maxiter; i++)
    {
        f(b, xn, args);
        matrix_mult_scalar(b, -1);
        Jacobian(J, xn, args);

        if(matrix_solve(J, b) != 0)
            break;

        matrix_add(xn, b, 1, NULL);

        if(eps > 0)
        {
            double norm;
            if(matrix_norm(b, 'F', &norm) < 0)
                break;

            if(norm < eps)
            {
                converged = 1;
                break;
            }
        }   
    }   

    matrix_free(J);
    matrix_free(b);

    if(converged)
        return i;
    else
        return -i;
}
