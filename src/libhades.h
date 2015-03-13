#ifndef LAPACK_HELPER_H
#define LAPACK_HELPER_H

    #include <complex.h>
    #include <stdio.h>

    #define LIBHADES_ERROR_OOM 254
    
    #ifndef M_PI
        #define M_PI 3.14159265358979323846
    #endif

    #ifndef M_E
        #define M_E 2.718281828459045
    #endif

    typedef double complex complex_t;
    typedef double complex complex_t;
    #define CREAL creal
    #define CIMAG cimag
    #define CABS cabs
    #define CONJ conj
    #define CPLX(a,b) (a+I*b)

    #define kronecker(i,j) (i==j)

    #include <libhades_clapack.h>
    #include <libhades_arpack.h>

    #define matrix_get(matrix,m,n)   (((matrix)->M)[(n)*((matrix)->rows)+(m)])
    #define matrix_set(matrix,m,n,v) (((matrix)->M)[(n)*((matrix)->rows)+(m)]=(v))

    #define pow_2(x) ((x)*(x))
    #define pow_3(x) ((x)*(x)*(x))

    /**
     * @brief implement MAX macro
     */
    #ifndef MAX
        #define MAX(a,b) \
          ({ __typeof__ (a) _a = (a); \
              __typeof__ (b) _b = (b); \
            _a > _b ? _a : _b; })
    #endif

    /**
     * @brief implement MIN macro
     */
    #ifndef MIN
        #define MIN(a,b) \
          ({ __typeof__ (a) _a = (a); \
              __typeof__ (b) _b = (b); \
            _a < _b ? _a : _b; })
    #endif

    #define LINSPACE(START, STOP, N, i) ((STOP-START)*((double)i/(N-1)))

    /**
     * @brief struct for real matrices
     */
    typedef struct {
        int rows;    /**< rows of matrix */
        int columns; /**< columns of matrix */
        int min;     /**< MIN(rows,columns) */
        int size;    /**< size = rows*columns */
        int type;    /**< type of matrix */
        int view;    /**< matrix is a view */
        double *M;   /**< data in column major order */
    } matrix_t;

    /**
     * @brief struct for complex matrices
     */
    typedef struct {
        int rows;       /**< rows of matrix */
        int columns;    /**< columns of matrix */
        int min;        /**< MIN(rows,columns) */
        int size;       /**< size = rows*columns */
        int type;       /**< type of matrix */
        int view;    /**< matrix is a view */
        complex_t *M;   /**< data in column major order */
    } matrix_complex_t;


    /* prototypes */
    int argmin   (double list[], int size);
    int argmax   (double list[], int size);
    int argabsmin(double list[], int size);

    void matrix_set_alloc(void *(*_malloc_cb)(size_t), void  (*_free_cb)(void *));

    matrix_t         *matrix_alloc        (int rows, int columns);
    matrix_complex_t *matrix_complex_alloc(int rows, int columns);

    matrix_t         *matrix_zeros        (int rows, int columns, matrix_t *A);
    matrix_complex_t *matrix_complex_zeros(int rows, int columns, matrix_complex_t *A);

    void matrix_setall        (matrix_t         *A, double x);
    void matrix_complex_setall(matrix_complex_t *A, complex_t x);

    matrix_t         *matrix_eye        (int dim, matrix_t *A);
    matrix_complex_t *matrix_complex_eye(int dim, matrix_complex_t *A);

    void matrix_free        (matrix_t         *matrix);
    void matrix_complex_free(matrix_complex_t *matrix);

    void matrix_fprintf        (FILE *stream, matrix_t         *A, const char *format, const char *sep, const char *sep_line);
    void matrix_complex_fprintf(FILE *stream, matrix_complex_t *A, const char *format, const char *sep, const char *sep_line);

    double matrix_trace_AB(matrix_t *A, matrix_t *B);
    complex_t matrix_trace_complex_AB(matrix_complex_t *A, matrix_complex_t *B);
    double matrix_trace_complex_AB_real(matrix_complex_t *A, matrix_complex_t *B);

    matrix_t         *matrix_kron        (matrix_t         *A, matrix_t         *B, matrix_t         *C);
    matrix_complex_t *matrix_complex_kron(matrix_complex_t *A, matrix_complex_t *B, matrix_complex_t *C);

    int matrix_add             (matrix_t         *A, matrix_t         *B, double    alpha, matrix_t         *C);
    int matrix_complex_add     (matrix_complex_t *A, matrix_complex_t *B, complex_t alpha, matrix_complex_t *C);
    int matrix_complex_add_real(matrix_complex_t *A, matrix_t         *B, complex_t alpha, matrix_complex_t *C);

    matrix_t         *matrix_copy        (matrix_t         *A);
    matrix_complex_t *matrix_complex_copy(matrix_complex_t *A);

    matrix_complex_t *matrix_tocomplex(matrix_t *A, matrix_complex_t *C);

    void matrix_transpose        (matrix_t         *A);
    void matrix_complex_transpose(matrix_complex_t *A);

    matrix_complex_t *matrix_mult_complex(matrix_t *A, complex_t alpha);

    void matrix_mult_scalar        (matrix_t         *A, double    alpha);
    void matrix_complex_mult_scalar(matrix_complex_t *A, complex_t alpha);

    matrix_complex_t *matrix_mult_complex_scalar(matrix_t *A, complex_t alpha, matrix_complex_t *C);

    matrix_t         *matrix_mult        (matrix_t *A,         matrix_t *B,         double    alpha, matrix_t *C);
    matrix_complex_t *matrix_complex_mult(matrix_complex_t *A, matrix_complex_t *B, complex_t alpha, matrix_complex_t *C);

    double    matrix_trace        (matrix_t *A);
    complex_t matrix_complex_trace(matrix_complex_t *A);

    matrix_t         *matrix_diag        (matrix_t *v);
    matrix_complex_t *matrix_complex_diag(matrix_complex_t *v);

    int eig_sym(matrix_t *A, char *JOBZ, char *UPLO, matrix_t *w);
    int eig_herm(matrix_complex_t *A, char *JOBZ, char *UPLO, matrix_t *w);
    int eig_complex_generic(matrix_complex_t *A, matrix_complex_t *w, matrix_complex_t *vl, matrix_complex_t *vr);

    int matrix_norm        (matrix_t         *A, char norm_type, double *norm);
    int matrix_complex_norm(matrix_complex_t *A, char norm_type, double *norm);

    matrix_complex_t *matrix_complex_exp_taylor(matrix_complex_t *A, int order);

    int matrix_lu_decomposition        (matrix_t         *A, int ipiv[]);
    int matrix_complex_lu_decomposition(matrix_complex_t *A, int ipiv[]);

    int matrix_invert        (matrix_t         *A);
    int matrix_complex_invert(matrix_complex_t *A);

    int rk4        (void (*f)(double t, matrix_t         *y, matrix_t         *ft, void *args), matrix_t         *yn, double t0, double t, double h, void *args);
    int rk4_complex(void (*f)(double t, matrix_complex_t *y, matrix_complex_t *ft, void *args), matrix_complex_t *yn, double t0, double t, double h, void *args);

    double vector_dot(matrix_t *x, matrix_t *y);
    complex_t vector_complex_dot(matrix_complex_t *x, matrix_complex_t *y);

    matrix_t         *matrix_get_column        (matrix_t         *A, int i);
    matrix_complex_t *matrix_complex_get_column(matrix_complex_t *A, int i);

    int sparse_complex_eig(int nx, int nev, char *which, void (*Av)(int nx, complex_t *in, complex_t *out, void *data), complex_t *d, int mxiter, double tol, void *data);

    int matrix_solve(matrix_t *A, matrix_t *b);

    int newton_mdim(void (*f)(matrix_t *, matrix_t *, void *), void (*Jacobian)(matrix_t *, matrix_t *, void *), matrix_t *xn, double eps, int maxiter, void *args);

#endif
