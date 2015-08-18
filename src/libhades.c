/**
 * @file   libhades.c
 * @author Michael Hartmann <michael.hartmann@physik.uni-augsburg.de>
 * @date   August, 2015
 * @brief  library to access low-level LAPACK functions
 */

#include <cblas.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>

#include <libhades.h>
#include <libhades/parse_npy_dict.h>

/** \defgroup misc miscellaneous functions
 *  @{
 */

/** @brief malloc wrapper
 *
 * This function uses malloc to allocate size bytes of memory and returns a
 * pointer to the memory. If an error occures, an error message is printed to
 * stderr and the program is aborted.
 *
 * @param [in] size amount of memory to allocate
 * @retval ptr pointer to the allocated memory
 */
void *xmalloc(size_t size)
{
    void *ptr = malloc(size);
    if(ptr == NULL)
    {   
        int err = errno;
        fprintf(stderr, "malloc can't allocate %zu bytes of memory: %s (%d)\n", size, strerror(err), err);
        abort();
    }   
    return ptr;
}

/** @brief free wrapper
 *
 * This function frees the memory allocated by \ref xmalloc or \ref xrealloc
 * (or malloc and realloc). If the pointer given is NULL, an error is printed
 * to stderr and the program is aborted.
 *
 * @param [in] ptr pointer to the memory that should be freed
 */
void xfree(void *ptr)
{
    if(ptr == NULL)
    {
        fprintf(stderr, "Trying to free NULL pointer\n");
        abort();
    }

    free(ptr);
}

/** @brief realloc wrapper
 *
 * This function changes the size of the memory block pointed to by ptr to size
 * bytes. If ptr is NULL, this function behaves like \ref xmalloc. If an error
 * occures, an error is printed to stderr and the program is aboprted.
 *
 * @param [in] ptr pointer to the memory block
 * @param [in] size size of the memory block
 * @retval ptr_new pointer to the new memory block
 */
void *xrealloc(void *ptr, size_t size)
{
    void *ptr_new = realloc(ptr, size);

    if(ptr_new == NULL)
    {
        int err = errno;
        fprintf(stderr, "realloc can't allocate %zu bytes of memory: %s (%d)\n", size, strerror(err), err);
        abort();
    }

    return ptr_new;
}

static void *(*malloc_cb)(size_t)          = &xmalloc;
//static void *(*realloc_cb)(void *, size_t) = &xrealloc;
static void  (*free_cb)(void *)            = &xfree;


/** macro to create functions argmin, argmax, argabsmin, argabsmax */
#define ARGXXX(FUNCTION_NAME, FUNCTION, RELATION) \
int FUNCTION_NAME(double list[], int size) \
{ \
    int index = 0; \
    for(int i = 0; i < size; i++) \
        if(FUNCTION(list[i]) RELATION FUNCTION(list[index])) \
            index = i; \
    return index; \
}

/** @brief Return index of smallest element in list
 *
 * @param [in] list
 * @param [in] size elements in list
 *
 * @retval index
 */
ARGXXX(argmin, +, <)

/** @brief Return index of largest element in list
 *
 * @param [in] list
 * @param [in] size elements in list
 *
 * @retval index
 */
ARGXXX(argmax, +, >)

/** @brief Return index of element with smallest absolute value in list
 *
 * @param [in] list
 * @param [in] size elements in list
 *
 * @retval index
 */
ARGXXX(argabsmin, fabs, <)
/** @brief Return index of element with largest absolute value in list
 *
 * @param [in] list
 * @param [in] size elements in list
 *
 * @retval index
 */
ARGXXX(argabsmax, fabs, >)

/** @}*/



/** \defgroup create Creating, printing and freeing matrices
 *  @{
 */

#define MATRIX_SWAP(FUNCTION_NAME, MATRIX_TYPE, TYPE) \
void FUNCTION_NAME(MATRIX_TYPE *A, MATRIX_TYPE *B) \
{ \
    int temp; \
    TYPE *ptr = A->M; \
    A->M = B->M; \
    B->M = ptr; \
    \
    temp = A->rows; \
    A->rows = B->rows; \
    B->rows = temp; \
    \
    temp = A->columns; \
    A->columns = B->columns; \
    B->columns = temp; \
}

/** @brief Swap matrices A and B
 *
 * This function swaps the matrices A and B. The former content of A will be
 * the content of B and vice versa. No data is copied or moved, but the
 * pointers are swapped.
 *
 * @param [in,out] A matrix A
 * @param [in,out] B matrix B
 */
MATRIX_SWAP(matrix_swap, matrix_t, double);

/** @brief Swap matrices A and B
 *
 * See \ref matrix_swap.
 *
 * @param [in,out] A matrix A
 * @param [in,out] B matrix B
 */
MATRIX_SWAP(matrix_complex_swap, matrix_complex_t, complex_t);

/** macro to create a diagnal matrix out of a vector */
#define MATRIX_DIAG(FUNCTION_NAME, MATRIX_TYPE, ZEROS) \
MATRIX_TYPE *FUNCTION_NAME(MATRIX_TYPE *v) \
{ \
    int dim = v->size; \
    MATRIX_TYPE *A = ZEROS(dim, dim, NULL); \
    if(A == NULL) \
        return NULL; \
\
    for(int i = 0; i < dim; i++) \
        matrix_set(A, i,i, v->M[i]); \
\
    return A; \
}

/** @brief Construct a real diagonal matrix from a row or columns vector
 *
 * @param [in] v row or column vector
 *
 * @retval A A = diag(v) if successfull, NULL otherwise
 */
MATRIX_DIAG(matrix_diag, matrix_t, matrix_zeros)

/** @brief Construct a complex diagonal matrix from a row or columns vector
 *
 * @param [in] v row or column vector
 *
 * @retval A A = diag(v) if successfull, NULL otherwise
 */
MATRIX_DIAG(matrix_complex_diag, matrix_complex_t, matrix_complex_zeros)


/** @brief Set functions to allocate and free memory
 *
 * By default wrappers to malloc and free from <stdlib.h> are used to allocate
 * and free memory. If allocation of memory fails or a NULL pointer is freed,
 * the program will terminate.
 *
 * @param [in] _malloc_cb callback to a malloc-alike function
 * @param [in] _free_cb callback to a free-alike function
 */
void matrix_set_alloc(void *(*_malloc_cb)(size_t), void  (*_free_cb)(void *))
{
    malloc_cb = _malloc_cb;
    free_cb   = _free_cb;
}


/** macro for copying matrices. */
#define MATRIX_COPY(FUNCTION_NAME, MTYPE, TYPE, ALLOC) \
MTYPE *FUNCTION_NAME(MTYPE *A, MTYPE *C) \
{ \
    if(C == NULL) { \
        C = ALLOC(A->rows, A->columns); \
        if(C == NULL) \
            return NULL; \
    } \
\
    C->rows    = A->rows; \
    C->columns = A->columns; \
    C->min     = A->min; \
    C->size    = A->size; \
    C->type    = A->type; \
    C->view    = 0; \
    memcpy(C->M, A->M, C->size*sizeof(TYPE)); \
    return C; \
}



/** @brief Copy real matrix A
 *
 * Copy matrix A into C. If C is NULL, space for the matrix C will be
 * allocated.
 *
 * @param [in]     A real matrix
 * @param [in,out] C real matrix
 *
 * @retval C copy of A
 */
MATRIX_COPY(matrix_copy, matrix_t, double, matrix_alloc)


/** @brief Copy complex matrix A
 *
 * Copy matrix A into C. If C is NULL, space for the matrix C will be
 * allocated.
 *
 * @param [in]     A complex matrix
 * @param [in,out] C complex matrix
 *
 * @retval C copy of A
 */
MATRIX_COPY(matrix_complex_copy, matrix_complex_t, complex_t, matrix_complex_alloc)


/** @brief Copy a real matrix A to a complex matrix C
 *
 * Copy the matrix A to a complex matrix C. The matrix elements of A and C will
 * be identical. If C is NULL, memory for the matrix will be allocated.
 *
 * @param [in] A real matrix
 *
 * @retval C copy of A if successfull, NULL otherwise
 */
matrix_complex_t *matrix_tocomplex(matrix_t *A, matrix_complex_t *C)
{
    if(C == NULL)
    {
        C = matrix_complex_alloc(A->rows, A->columns);
        if(C == NULL)
            return NULL;
    }

    for(int i = 0; i < C->size; i++)
        C->M[i] = A->M[i];

    return C;
}

/** @brief Print real matrix M to stream
 *
 * Print the matrix A to the stream given by stream, e.g. stdout or stderr.
 * The format is given by the format string format, the separator of two
 * columns is sep, the separator between lines is given by sep_line. Both sep
 * and sep_line may be NULL.
 *
 * @param [in] stream output stream
 * @param [in] A real matrix
 * @param [in] format output format, e.g. "%lf" or "%g"
 * @param [in] sep separator between columns, e.g. "\t"
 * @param [in] sep_line separator between lines, e.g. "\n"
 */
void matrix_fprintf(FILE *stream, matrix_t *A, const char *format, const char *sep, const char *sep_line)
{
    const int rows = A->rows, columns = A->columns;

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            fprintf(stream, format, matrix_get(A, i,j));
            if(sep != NULL)
                fputs(sep, stream);
        }
        if(sep_line != NULL)
            fputs(sep_line, stream);
    }
}

/** @brief Print complex matrix A to stream
 *
 * See matrix_fprintf.
 *
 * @param [in] stream output stream
 * @param [in] A complex matrix
 * @param [in] format output format, e.g. "%+lf%+lfi"
 * @param [in] sep separator between columns, e.g. "\t"
 * @param [in] sep_line separator between lines, e.g. "\n"
 */
void matrix_complex_fprintf(FILE *stream, matrix_complex_t *A, const char *format, const char *sep, const char *sep_line)
{
    const int rows = A->rows, columns = A->columns;

    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
        {
            const complex_t c = matrix_get(A, i,j);
            fprintf(stream, format, CREAL(c), CIMAG(c));
            if(sep != NULL)
                fputs(sep, stream);
        }
        if(sep_line != NULL)
            fputs(sep_line, stream);
    }
}

/** macro for matrix allocations. */
#define MATRIX_ALLOC(FUNCTION_NAME, MTYPE, TYPE) \
MTYPE *FUNCTION_NAME(int rows, int columns) \
{ \
    MTYPE *A = malloc_cb(sizeof(MTYPE)); \
    if(A == NULL) \
        return NULL; \
\
    A->rows    = rows; \
    A->columns = columns; \
    A->min     = MIN(rows, columns); \
    A->size    = rows*columns; \
    A->type    = 0; \
    A->view    = 0; \
    A->M       = malloc_cb(rows*columns*sizeof(TYPE)); \
    if(A->M == NULL) \
    { \
        free_cb(A); \
        return NULL; \
    } \
\
    return A; \
}

/** @brief Allocate memory for real matrix of m rows and n columns
 *
 * This function will allocate and return a matrix of m lines and n columns.
 * The function given by matrix_set_alloc will be used to allocate memory; by
 * default this is malloc from <stdlib.h>.
 *
 * The matrix elements will be undefined. To allocate a matrix initialized with
 * zeros, see matrix_zeros. To create a unity matrix, see matrix_eye.
 *
 * @param [in] rows rows of matrix M
 * @param [in] columns columns of matrix M
 *
 * @retval A real matrix if successful, NULL otherwise
 */
MATRIX_ALLOC(matrix_alloc, matrix_t, double)

/** @brief Allocate memory for complex matrix of m rows and n columns
 *
 * See matrix_alloc.
 *
 * @param [in] rows rows of matrix M
 * @param [in] columns columns of matrix M
 *
 * @retval A complex matrix if successful, otherwise NULL
 */
MATRIX_ALLOC(matrix_complex_alloc, matrix_complex_t, complex_t)

/** macro to create zero matrix */
#define MATRIX_ZEROS(FUNCTION_NAME, MTYPE, TYPE, ALLOC, SETALL) \
MTYPE *FUNCTION_NAME(int rows, int columns, MTYPE *A) \
{ \
    if(A == NULL) \
    { \
        A = ALLOC(rows,columns); \
        if(A == NULL) \
            return NULL; \
    } \
\
    SETALL(A, 0); \
\
    return A; \
}

/** @brief Generate a real zero matrix of m lines and n columns
 *
 * Set every element of matrix A to 0. If A is NULL, the matrix will be
 * created. In this case, you have to free the matrix yourself.
 *
 * @param [in] rows rows of matrix M
 * @param [in] columns columns of matrix M
 * @param [in,out] A: matrix
 *
 * @retval A real matrix 0 if successful, otherwise NULL
 */
MATRIX_ZEROS(matrix_zeros, matrix_t, double, matrix_alloc, matrix_setall)

/** @brief Generate a complex zero matrix of m lines and n columns
 *
 * See matrix_zeros.
 *
 * @param [in] rows rows of matrix M
 * @param [in] columns columns of matrix M
 * @param [in,out] A: matrix
 *
 * @retval A complex matrix 0 if successful, otherwise NULL
 */
MATRIX_ZEROS(matrix_complex_zeros, matrix_complex_t, complex_t, matrix_complex_alloc, matrix_complex_setall)


/** macro to create matrix x*Id */
#define MATRIX_SETALL(FUNCTION_NAME, MATRIX_TYPE, TYPE) \
void FUNCTION_NAME(MATRIX_TYPE *A, TYPE x) \
{ \
    const int size = A->size; \
    TYPE *M = A->M; \
    for(int i = 0; i < size; i++) \
        *M++ = x; \
}

/** @brief Set matrix elements of A to x
 *
 * @param [in,out] A real matrix
 * @param [in] x real number
 */
MATRIX_SETALL(matrix_setall, matrix_t, double)

/** @brief Set matrix elements of A to x
 *
 * @param [in,out] A complex matrix
 * @param [in] x complex number
 */
MATRIX_SETALL(matrix_complex_setall, matrix_complex_t, complex_t)

/** macro to create unity matrix */
#define MATRIX_EYE(FUNCTION_NAME, MTYPE, TYPE, ALLOC, SETALL) \
MTYPE *FUNCTION_NAME(int dim, MTYPE *A) \
{ \
    int min; \
    TYPE *M; \
\
    if(A == NULL) \
    { \
        A = ALLOC(dim,dim); \
        if(A == NULL) \
            return NULL; \
    } \
 \
    min = A->min; \
    M = A->M; \
    SETALL(A,0); \
    for(int i = 0; i < min; i++) \
        *(M+i*(min+1)) = 1; \
 \
    return A; \
}


/** @brief Create real identity matrix
 *
 * If A == NULL, a identity matrix of dimension dim times dim is created and
 * returned.
 * If A != NULL, the matrix A is set to the identity matrix and the parameter
 * dim is ignored. More specific, for a general (e.g. not square) matrix A, the
 * matrix elements are set to A_ij = Delta_ij.
 *
 * @param [in]     dim dimension of identity matrix (ignored for A == NULL)
 * @param [in,out] A real matrix
 *
 * @retval A identity matrix if successful, NULL otherwise
 */
MATRIX_EYE(matrix_eye, matrix_t, double, matrix_alloc, matrix_setall)

/** @brief Create complex identity matrix
 *
 * See matrix_eye.
 *
 * @param [in]     dim dimension of identity matrix (ignored for A == NULL)
 * @param [in,out] A complex matrix
 *
 * @retval A identity matrix if successful, NULL otherwise
 */
MATRIX_EYE(matrix_complex_eye, matrix_complex_t, complex_t, matrix_complex_alloc, matrix_complex_setall)


/** macro to free matrices */
#define MATRIX_FREE(FUNCTION_NAME, MTYPE) \
void FUNCTION_NAME(MTYPE *A) \
{ \
    if(A != NULL) \
    { \
        if(!A->view && A->M != NULL) \
        { \
            free_cb(A->M); \
            A->M = NULL; \
        } \
        free_cb(A); \
    } \
}

/** @brief Free real matrix
 *
 * This function will free the memory allocated for matrix M. If M is NULL this
 * function will do nothing.
 *
 * @param [in,out] A matrix to free
 */
MATRIX_FREE(matrix_free, matrix_t)

/** @brief Free complex matrix
 *
 * See matrix_free.
 *
 * @param [in,out] A matrix to free
 */
MATRIX_FREE(matrix_complex_free, matrix_complex_t)

/** @}*/


/** \defgroup trdet trace and determinant
 *  @{
 */

/** macro to calculate trace of matrix */
#define MATRIX_TRACE(FUNCTION_NAME, MATRIX_TYPE, TYPE) \
TYPE FUNCTION_NAME(MATRIX_TYPE *A) \
{ \
    const int min = A->min; \
    const TYPE *M = A->M; \
    TYPE trace = 0; \
    for(int i = 0; i < min; i++) \
        trace += *(M+i*(1+min)); \
\
    return trace; \
}

/** @brief Calculate Tr(A) of real matrix A
 *
 * @param [in] A complex matrix
 *
 * @retval x with x=Tr(A)
 */
MATRIX_TRACE(matrix_trace, matrix_t, double)

/** @brief Calculate Tr(A) of complex matrix A
 *
 * @param [in] A complex matrix
 *
 * @retval z with z=Tr(A)
 */
MATRIX_TRACE(matrix_complex_trace, matrix_complex_t, complex_t)

/** @brief Calculate Tr(A*B)
 *
 * This will calculate the trace of A*B: Tr(A*B). Matrix A and B must be square
 * matrices of dimension dim, dim.
 *
 * A and B must not point to the same matrix or the behaviour will be
 * undefined!
 *
 * @param [in] A real matrix
 * @param [in] B real matrix
 *
 * @retval x with x=Tr(A*B)
 */
double matrix_trace_AB(matrix_t *A, matrix_t *B)
{
    const int dim = A->rows;
    double sum = 0;
    double *M1 = A->M;
    double *M2 = B->M;

    for(int i = 0; i < dim; i++)
        sum += cblas_ddot(dim, M1+i, dim, M2+i*dim, 1);

    return sum;
}

/** @brief Calculate Tr(A*B) for A,B complex
 *
 * See matrix_trace_AB.
 *
 * @param [in] A complex matrix
 * @param [in] B complex matrix
 *
 * @retval z with z=Tr(A*B)
 */
complex_t matrix_trace_complex_AB(matrix_complex_t *A, matrix_complex_t *B)
{
    const int dim = A->rows;
    complex_t sum = 0;

    for(int i = 0; i < dim; i++)
        for(int k = 0; k < dim; k++)
            sum += matrix_get(A, i,k)*matrix_get(B, k,i);

    return sum;
}

/** @brief Calculate Re(Tr(A*B)) for A,B complex
 *
 * See matrix_trace_AB.
 *
 * @param [in] A complex matrix
 * @param [in] B complex matrix
 *
 * @retval x with x=Re(Tr(A*B))
 */
double matrix_trace_complex_AB_real(matrix_complex_t *A, matrix_complex_t *B)
{
    const int dim = A->rows;
    double sum = 0;
    const complex_t *M1 = A->M;
    const complex_t *M2 = B->M;

    for(int i = 0; i < dim; i++)
        for(int k = 0; k < dim; k++)
            sum += CREAL(M1[i*dim+k]*M2[k*dim+i]);

    return sum;
}

/** @}*/


/** \defgroup kron Kronecker product
 *  @{
 */

/** macro to create Kronecker product */
#define MATRIX_KRON(FUNCTION_NAME, MTYPE, TYPE, ALLOC, SETALL) \
MTYPE *FUNCTION_NAME(MTYPE *A, MTYPE *B, MTYPE *C) \
{ \
    const int Am = A->rows, An = A->columns; \
    const int Bm = B->rows, Bn = B->columns; \
    if(C == NULL) \
    { \
        C = ALLOC(Am*Bm, An*Bn); \
        if(C == NULL) \
            return NULL; \
    } \
    SETALL(C, 0); \
\
    for(int m = 0; m < Am; m++) \
        for(int n = 0; n < An; n++) \
        { \
            const TYPE c = matrix_get(A,m,n); \
            if(c != 0) \
            { \
                for(int im = 0; im < Bm; im++) \
                    for(int in = 0; in < Bn; in++) \
                        matrix_set(C, m*Bm+im, n*Bn+in, c*matrix_get(B,im,in)); \
            } \
        } \
\
    return C; \
}


/** @brief Calculate Kronecker product of real matrices A and B
 *
 * Matrix A has dimension Am,An, matrix B has dimension Bm,Bn.
 *
 * If C is not NULL, the Kronecker product will be stored in C. C must have
 * dimension Am+Bm,An+Bn. If C is NULL, memory for the matrix will be allocated
 * and the matrix will be returned. You have to free the memory for the
 * returned matrix yourself.
 *
 * @param [in] A real matrix
 * @param [in] B real matrix
 * @param [in,out] C Kronecker product
 *
 * @retval C Kronecker product of A and B, NULL if no memory could be allocated
 */
MATRIX_KRON(matrix_kron, matrix_t, double, matrix_alloc, matrix_setall)


/** @brief Calculate Kronecker product of complex matrices A and B
 *
 * See matrix_kron.
 *
 * @param [in] A complex matrix
 * @param [in] B complex matrix
 * @param [in,out] C Kronecker product
 *
 * @retval C Kronecker product of A and B
 */
MATRIX_KRON(matrix_complex_kron, matrix_complex_t, complex_t, matrix_complex_alloc, matrix_complex_setall)

/** @}*/


/** \defgroup addsubmult Add, subtract and multiply matrices
 *  @{
 */


/** macro to multiply matrix with scalar factor */
#define MATRIX_MULT_SCALAR(FUNCTION_NAME, MTYPE, TYPE) \
void FUNCTION_NAME(MTYPE *A, TYPE alpha) \
{ \
    const int max = A->rows*A->columns; \
    TYPE *M = A->M; \
    for(int i = 0; i < max; i++) \
        M[i] *= alpha; \
}

/** @brief Multiply real matrix with a real scalar
 *
 * alpha*A -> A
 *
 * @param [in,out] A real matrix
 * @param [in] alpha real scalar
 * @retval A, A=alpha*A
 */
MATRIX_MULT_SCALAR(matrix_mult_scalar, matrix_t, double)

/** @brief Multiply complex matrix with a complex scalar
 *
 * alpha*A -> A
 *
 * @param [in,out] A complex matrix
 * @param [in] alpha complex scalar
 * @retval A, A=alpha*A
 */
MATRIX_MULT_SCALAR(matrix_complex_mult_scalar, matrix_complex_t, complex_t)

/** @brief Multiply real matrix with a complex scalar
 *
 * alpha*A -> C
 *
 * If C is NULL, memory for the matrix C will be allocated.
 *
 * @param [in] A real matrix
 * @param [in] alpha complex scalar
 * @param [in,out] C complex matrix
 *
 * @retval C, C = alpha*A
 */
matrix_complex_t *matrix_mult_complex_scalar(matrix_t *A, complex_t alpha, matrix_complex_t *C)
{
    const int rows    = A->rows;
    const int columns = A->columns;
    const double *AM     = A->M;
    complex_t *CM;
    if(C == NULL)
    {
        C = matrix_complex_alloc(rows, columns);
        if(C == NULL)
            return NULL;
    }
    CM = C->M;

    for(int i = 0; i < rows*columns; i++)
        CM[i] = alpha*AM[i];

    return C;
}


/* compute alpha*A*B */
#define MATRIX_MULT(FUNCTION_NAME, MATRIX_TYPE, TYPE, ALLOC, BLAS_xGEMM, PTR) \
MATRIX_TYPE *FUNCTION_NAME(MATRIX_TYPE *A, MATRIX_TYPE *B, TYPE alpha, MATRIX_TYPE *C) \
{ \
    TYPE beta = 0; \
\
    if(A->columns != B->rows) \
        return NULL; \
\
    if(C == NULL) \
    { \
        C = ALLOC(A->rows, B->columns); \
        if(C == NULL) \
            return NULL; \
    } \
\
    /* xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) */ \
    BLAS_xGEMM(CblasColMajor, /* column order */ \
               CblasNoTrans,  /* don't transpose/conjugate A */ \
               CblasNoTrans,  /* don't transpose/conjugate B */ \
               A->rows,       /* M: rows of A and C */ \
               B->columns,    /* N: columns of B and C */ \
               A->columns,    /* K: columns of A and rows of B */ \
               PTR alpha,     /* alpha: scalar */ \
               A->M,          /* A: matrix A */ \
               A->rows,       /* LDA: leading dimension of A (columns) */ \
               B->M,          /* B: matrix B */ \
               A->columns,    /* LDB: leading dimension of B (columns) */ \
               PTR beta,      /* beta: scalar */ \
               C->M,          /* C: matrix C */ \
               C->rows        /* LDC: leading dimension of C (columns) */ \
    ); \
\
    return C; \
}

/** @brief Multiply real matrices
 *
 * alpha*A*B -> C
 *
 * If C is NULL, memory for the matrix C will be allocated.
 *
 * @param [in] A real matrix
 * @param [in] B real matrix
 * @param [in] alpha real scalar
 * @param [in,out] C real matrix
 *
 * @retval C, C = alpha*A*B
 */
MATRIX_MULT(matrix_mult, matrix_t, double, matrix_alloc, cblas_dgemm, +)

/** @brief Multiply complex matrices
 *
 * alpha*A*B -> C
 *
 * If C is NULL, memory for the matrix C will be allocated.
 *
 * @param [in] A complex matrix
 * @param [in] B complex matrix
 * @param [in] alpha complex scalar
 * @param [in,out] C complex matrix
 *
 * @retval C, C = alpha*A*B
 */
MATRIX_MULT(matrix_complex_mult, matrix_complex_t, complex_t, matrix_complex_alloc, cblas_zgemm, &)

/** macro to add two matrices */
#define MATRIX_ADD(FUNCTION_NAME, TYPE1, MTYPE1, TYPE2, MTYPE2) \
int FUNCTION_NAME(MTYPE1 *A, MTYPE2 *B, TYPE1 alpha, MTYPE1 *C) \
{ \
    const int max = A->size; \
    TYPE1 *M3; \
    TYPE1 *M1 = A->M; \
    TYPE2 *M2 = B->M; \
\
    if(C == NULL) \
        M3 = A->M; \
    else \
        M3 = C->M; \
\
    if(A->rows != B->rows || A->columns != B->columns) \
        return LIBHADES_ERROR_SHAPE; \
\
    for(int i = 0; i < max; i++) \
        M3[i] = M1[i] + (alpha*M2[i]); \
\
    return 0; \
}

/** @brief Add real matrices A and B
 *
 * Calculate A+alpha*B -> C.
 *
 * The result will be stored in C. If C is NULL, the result is stored in A.
 *
 * @param [in,out] A real matrix
 * @param [in] B real matrix
 * @param [in] alpha real scalar
 * @param [in,out] C real matrix or NULL
 *
 * @retval 0 if successfull
 * @retval LIBHADES_ERROR_SHAPE if matrices have wrong shape
 */
MATRIX_ADD(matrix_add, double, matrix_t, double, matrix_t)

/** @brief Add complex matrices A and B
 *
 * Calculate A+alpha*B -> C.
 *
 * The result will be stored in C. If C is NULL, the result is stored in A.
 *
 * @param [in,out] A complex matrix
 * @param [in] B complex matrix
 * @param [in] alpha complex number
 * @param [in,out] C complex matrix or NULL
 *
 * @retval 0 if successfull
 * @retval LIBHADES_ERROR_SHAPE if matrices have wrong shape
 */
MATRIX_ADD(matrix_complex_add, complex_t, matrix_complex_t, complex_t, matrix_complex_t)

/** @brief Add complex matrix A and real matrix B
 *
 * Calculate A+alpha*B -> C.
 *
 * The result will be stored in C. If C is NULL, the result is stored in A.
 *
 * @param [in,out] A complex matrix
 * @param [in] B real matrix
 * @param [in] alpha complex scalar
 * @param [in,out] C complex matrix or NULL
 *
 * @retval 0 if successfull
 * @retval LIBHADES_ERROR_SHAPE if matrices have wrong shape
 */
MATRIX_ADD(matrix_complex_add_real, complex_t, matrix_complex_t, double, matrix_t)

/** @}*/


/** \defgroup transconj Transpose, conjugate
 *  @{
 */


#define MATRIX_TRANSPOSE(FUNCTION_NAME, MTYPE, TYPE) \
void FUNCTION_NAME(MTYPE *A) \
{ \
    const int rows    = A->rows; \
    const int columns = A->columns; \
    for(int im = 0; im < rows; im++) \
        for(int in = im+1; in < columns; in++) \
        { \
            TYPE temp = matrix_get(A, im, in); \
            matrix_set(A, im, in, matrix_get(A, in, im)); \
            matrix_set(A, in, im, temp); \
        } \
\
    A->rows    = columns; \
    A->columns = rows; \
}

/** @brief Transpose real matrix A
 *
 * The matrix will be transposed.
 *
 * @param [in,out] A real matrix
 *
 * @retval C with C=A^T
 */
MATRIX_TRANSPOSE(matrix_transpose, matrix_t, double)

/** @brief Transpose complex matrix A
 *
 * The matrix will be transposed.
 *
 * @param [in,out] A complex matrix
 *
 * @retval C with C=A^T
 */
MATRIX_TRANSPOSE(matrix_complex_transpose, matrix_complex_t, complex_t)

/** @}*/


/** \defgroup ev Eigenvalue problems
 *  @{
 */


/** @brief Compute eigenvalues and optionally eigenvectors of symmetric matrix A
 *
 * This function computes all eigenvalues and, optionally, eigenvectors of a
 * real symmetric matrix A.
 *
 * See dsyev.
 *
 * @param [in] A real matrix
 * @param [in] JOBZ 'N': only eigenvalues, 'V' eigenvalues and eigenvectors
 * @param [in] UPLO 'U': upper triangle part of A is stored; 'L': lower triangle part of A is stored
 * @param [in] w real matrix of dimension (dim,1) (i.e. a vector); the eigenvalues will be stored in w
 *
 * @retval 0 on success
 */
int eig_sym(matrix_t *A, char *JOBZ, char *UPLO, matrix_t *w)
{
    int info, lwork = -1, N = A->min;
    double workopt;
    double *work;

    dsyev_(JOBZ, UPLO, &N, A->M, &N, w->M, &workopt, &lwork, &info);
    if(info != 0)
        return info;

    lwork = workopt;
    work = malloc_cb(lwork*sizeof(double));
    if(work == NULL)
        return LIBHADES_ERROR_OOM;

    dsyev_(JOBZ, UPLO, &N, A->M, &N, w->M, work, &lwork, &info);

    free_cb(work);

    return info;
}

/** @brief Compute eigenvalues and optionally eigenvectors of Hermitian matrix A
 *
 * This function computes all eigenvalues and, optionally, eigenvectors of a
 * Hermitian symmetric matrix A.
 *
 * See zheev.
 *
 * @param [in] A complex matrix
 * @param [in] JOBZ 'N': only eigenvalues, 'V' eigenvalues and eigenvectors
 * @param [in] UPLO 'U': upper triangle part of A is stored; 'L': lower triangle part of A is stored
 * @param [in] w real matrix of dimension (dim,1) (i.e. a vector); the eigenvalues will be stored in w
 *
 * @retval 0 on success
 */
int eig_herm(matrix_complex_t *A, char *JOBZ, char *UPLO, matrix_t *w)
{
    int info, lwork = -1, N = A->min;
    complex_t workopt;
    complex_t *work = NULL;
    double *rwork;
    
    rwork = malloc_cb(MAX(1, 3*N-2)*sizeof(double));
    if(rwork == NULL)
        return LIBHADES_ERROR_OOM;

    zheev_(JOBZ, UPLO, &N, A->M, &N, w->M, &workopt, &lwork, rwork, &info);
    if(info != 0)
    {
        free_cb(rwork);
        return info;
    }

    lwork = CREAL(workopt);
    work = malloc_cb(lwork*sizeof(complex_t));
    if(work == NULL)
    {
        free_cb(rwork);
        return LIBHADES_ERROR_OOM;
    }

    zheev_(JOBZ, UPLO, &N, A->M, &N, w->M, work, &lwork, rwork, &info);

    free_cb(rwork);
    free_cb(work);

    return info;
}

/** @brief Compute eigenvalues and optionally eigenvectors of matrix A
 *
 * Compute for an N-by-N complex nonsymmetric matrix A, the eigen- values and,
 * optionally, right eigenvectors
 *
 * See zgeev.
 *
 * @param [in] A real matrix
 * @param [in] w list containing the eigenvalues of A
 * @param [in] vr if vr != NULL, right eigenvectors are computed and will be stored in vr; vr must be a complex matrix of dimension (dim,1) (i.e. a vector)
 * @param [in] vl if vl != NULL, lefft eigenvectors are computed and will be stored in vl; vl must be a complex matrix of dimension (dim,1) (i.e. a vector)
 *
 * @retval 0 on success
 */

/* Parameters */
int eig_complex_generic(matrix_complex_t *A, matrix_complex_t *w, matrix_complex_t *vl, matrix_complex_t *vr)
{
    int N = A->min;
    int lwork = -1;
    char *jobvr = "N";
    char *jobvl = "N";
    int info;
    complex_t *evr  = NULL;
    complex_t *evl  = NULL;
    complex_t *work = NULL;
    double *rwork   = NULL;
    complex_t wopt;

    /* if vr is not NULL, calculate right eigenvectors */
    if(vr != NULL)
    {
        jobvr = "V";
        evr   = vr->M;
    }
    /* if vl is not NULL, calculate left eigenvectors */
    if(vl != NULL)
    {
        jobvl = "V";
        evl   = vl->M;
    }

    /* get the optimal size for workspace work */
    zgeev_(jobvl, jobvr, &N, A->M, &N, w->M, evl, &N, evr, &N, &wopt, &lwork, rwork, &info);

    if(info != 0)
        return info;

    lwork = CREAL(wopt);

    rwork = malloc_cb(2*N*sizeof(double));
    work  = malloc_cb(lwork*sizeof(complex_t));

    if(rwork == NULL || work == NULL)
    {
        if(rwork != NULL)
            free(rwork);
        if(work != NULL)
            free(work);

        return LIBHADES_ERROR_OOM;
    }

    /* SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO ) */
    zgeev_(
        jobvl,      /* left eigenvectors of A are not computed */
        jobvr,      /* calculate/don't calculate right eigenvectors */
        &N,         /* order of matrix A */
        A->M,       /* matrix A */
        &N,         /* LDA - leading dimension of A */
        w->M,       /* eigenvalues */
        evl,        /* left eigenvectors */
        &N,         /* leading dimension of array vl */
        evr,        /* eigenvectors */
        &N,         /* leading dimension of the array VR */
        work,       /* COMPLEX*16 array, dimension (LWORK) */
        &lwork,     /* dimension of the array WORK; LWORK >= max(1,2*N) */
        rwork,      /* (workspace) DOUBLE PRECISION array, dimension (2*N) */
        &info       /* 0 == success */
    );

    free_cb(work);
    free_cb(rwork);

    return info;
}

/** @}*/




/** \defgroup exp Calculate matrix exponential
 *  @{
 */

/** macro to calculate matrix norm */
#define MATRIX_NORM(FUNCTION_NAME, MTYPE, LAPACK_FUNC) \
int FUNCTION_NAME(MTYPE *A, char norm_type, double *norm) \
{ \
    double *work = NULL; \
\
    if(norm_type == 'I') \
    { \
        work = malloc_cb(A->rows*sizeof(double)); \
        if(work == NULL) \
            return LIBHADES_ERROR_OOM; \
    } \
\
    *norm = LAPACK_FUNC(&norm_type, &A->rows, &A->columns, A->M, &A->rows, work); \
\
    if(work != NULL) \
        free_cb(work); \
\
    return 0; \
}

/** @brief Compute matrix norm for real matrix
 *
 * See dlange.
 *
 * Returns
 *      max(abs(A(i,j))), norm_type = 'M' or 'm'
 *      norm1(A),         norm_type = '1', 'O' or 'o'
 *      normI(A),         norm_type = 'I' or 'i'
 *      normF(A),         norm_type = 'F', 'f', 'E' or 'e'
 *
 * @param [in]  A real matrix
 * @param [in]  norm_type type of norm, e.g. 'F' for Frobenius norm
 * @param [out] norm value of norm
 *
 * @retval ret 0 if successfull, <0 otherwise
 */
MATRIX_NORM(matrix_norm, matrix_t, dlange_);

/** @brief Compute matrix norm for complex matrix
 *
 * See matrix_norm.
 *
 * @param [in]  A complex matrix
 * @param [in]  norm_type type of norm, e.g. 'F' for Frobenius norm
 * @param [out] norm value of norm
 *
 * @retval ret 0 if successfull, <0 otherwise
 */
MATRIX_NORM(matrix_complex_norm, matrix_complex_t, zlange_);

matrix_complex_t *matrix_complex_exp_taylor(matrix_complex_t *A, int order)
{
    const int rows = A->min;
    matrix_complex_t *C = matrix_complex_alloc(rows,rows);

    /* B = E+A/order */
    matrix_complex_t *B = matrix_complex_copy(A,NULL);
    matrix_complex_mult_scalar(B, 1./order);

    for(int i = 0; i < rows; i++)
        matrix_set(B, i,i, 1+matrix_get(B,i,i));

    for(int k = order-1; k > 0; k--)
    {
        matrix_complex_mult(A, B, 1./order, C);

        /* swap */
        {
            matrix_complex_t *temp;
            temp = C;
            C = B;
            B = temp;
        }

        for(int i = 0; i < rows; i++)
            matrix_set(B, i,i, 1+matrix_get(B,i,i));
    }

    matrix_complex_free(C);

    return B;
}

/** @}*/



/** \defgroup LA LU decomposition, inverting
 *  @{
 */

#define LU_DECOMPOSITION(FUNCTION_NAME, MATRIX_TYPE, XGETRF) \
int FUNCTION_NAME(MATRIX_TYPE *A, int ipiv[]) \
{ \
    int info; \
\
    XGETRF( \
        &A->rows,    /* M number of rows of A */ \
        &A->columns, /* N number of columns of A */ \
        A->M,        /* matrix A to be factored */ \
        &A->columns, /* LDA: leading dimension of A */ \
        ipiv,        /* pivot indices of dimension (min(M,N)) */ \
        &info \
    ); \
\
    return info; \
}

/** @brief Compute LU decomposition of real matrix A
 *
 * See dgetrf.
 *
 * The factorization has the form
 *    A = P * L * U
 * where P is a permutation matrix, L is lower triangular with unit
 * diagonal elements (lower trapezoidal if rows > columns), and U is upper
 * triangular (upper trapezoidal if m < n).
 *
 * @param [in,out] A real matrix
 * @param [out] ipiv pivot indices; array of dimension MIN(rows,columns)
 *
 * @retval INFO
 */
LU_DECOMPOSITION(matrix_lu_decomposition, matrix_t, dgetrf_)

/** @brief Compute LU decomposition of complex matrix A
 *
 * See matrix_lu_decomposition.
 *
 * @param [in,out] A complex matrix
 * @param [out] ipiv pivot indices; array of dimension MIN(rows,columns)
 *
 * @retval INFO
 */
LU_DECOMPOSITION(matrix_complex_lu_decomposition, matrix_complex_t, zgetrf_)

#define MATRIX_INVERT(FUNCTION_NAME, TYPE, MATRIX_TYPE, LU_DECOMPOSITION, XGETRI) \
int FUNCTION_NAME(MATRIX_TYPE *A) \
{ \
    int info, lwork, dim = A->min; \
    int *ipiv = NULL; \
    TYPE *work = NULL; \
    TYPE workopt; \
\
    ipiv = malloc_cb(dim*sizeof(int)); \
    if(ipiv == NULL) \
        return LIBHADES_ERROR_OOM; \
\
    info = LU_DECOMPOSITION(A, ipiv); \
    if(info != 0) \
        goto out; \
\
    lwork = -1; \
    XGETRI(&dim, A->M, &dim, ipiv, &workopt, &lwork, &info); \
    if(info != 0) \
        goto out; \
\
    lwork = (int)workopt; \
    work = malloc_cb(lwork*sizeof(TYPE)); \
    if(work == NULL) \
    { \
        info = LIBHADES_ERROR_OOM; \
        goto out; \
    } \
\
    XGETRI( \
        &dim,   /* order of matrix A */ \
        A->M,   /* factors L and U from LU decomposition */ \
        &dim,   /* LDA: leading dimension of A */ \
        ipiv,   /* pivot indices */ \
        work,   /* workspace of dimension LWORK */ \
        &lwork, /* length of work */ \
        &info \
    ); \
\
out: \
    free_cb(ipiv); \
    if(work != NULL) \
        free_cb(work); \
\
    return info; \
}

/** @brief Invert real matrix A
 *
 * The inverse of A is computed using the LU factorzation of A
 *
 * @param [in] A real matrix
 *
 * @retval INFO
 */
MATRIX_INVERT(matrix_invert, double, matrix_t, matrix_lu_decomposition, dgetri_)

/** @brief Invert complex matrix A
 *
 * See matrix_invert.
 *
 * @param [in] A complex matrix
 *
 * @retval INFO
 */
MATRIX_INVERT(matrix_complex_invert, complex_t, matrix_complex_t, matrix_complex_lu_decomposition, zgetri_)


/** @brief Solve system of linear equations
 *
 * Solve the system of linear equations:
 *      A*x = b
 *
 * @param [in,out] A matrix
 * @param [in,out] b vector/matrix
 *
 * @retval INFO
 */
int matrix_solve(matrix_t *A, matrix_t *b)
{
    int N = A->min;
    char trans = 'N';
    int nrhs = b->columns;
    int info = -1;
    int *ipiv = malloc_cb(N*sizeof(int));

    if(ipiv == NULL)
        return LIBHADES_ERROR_OOM;

    matrix_lu_decomposition(A, ipiv);
    dgetrs_(&trans, &N, &nrhs, A->M, &N, ipiv, b->M, &N, &info);

    free_cb(ipiv);

    return info;
}

/** @}*/

/*
static void _cblas_zaxpy(const int N, const double alpha, const void *X, const int incX, void *Y, const int incY)
{
    complex_t beta = alpha;
    cblas_zaxpy(N, &beta, X, incX, Y, incY);
}
*/

/** @brief Calculate dot product of vectors x and y
 *
 * Calculate dot product of first column of x,y. If x and y have different
 * rows, the minimum is used.
 *
 * @param [in] x vector
 * @param [in] y vector
 * @retval x*y
 */
double vector_dot(matrix_t *x, matrix_t *y)
{
    int incx = 1;
    int incy = 1;
    int N = MIN(x->rows, y->rows);
    
    return ddot_(&N, x->M, &incx, y->M, &incy);
}

/** @brief Calculate dot product of vectors x and y
 *
 * Calculate dot product of first column of x,y. If x and y have different
 * rows, the minimum is used.
 *
 * @param [in] x vector
 * @param [in] y vector
 * @retval x*y
 */
complex_t vector_complex_dot(matrix_complex_t *x, matrix_complex_t *y)
{
    /*
    int incx = 1;
    int incy = 1;
    int N = MIN(x->rows, y->rows);
    complex_t z = 0;
    
    zdotc_(&z, &N, x->M, &incx, y->M, &incy);
    return z;
    */

    int N = MIN(x->rows, y->rows);
    complex_t z = 0;
    for(int i = 0; i < N; i++)
        z += matrix_get(x,i,0)*matrix_get(y,i,0);
    return z;
}

#define MATRIX_GET_COLUMN(FUNCTION_NAME, TYPE, MATRIX_TYPE) \
MATRIX_TYPE *FUNCTION_NAME(MATRIX_TYPE *A, int i) \
{ \
    const int rows = A->rows; \
    MATRIX_TYPE *v = malloc_cb(sizeof(MATRIX_TYPE)); \
    if(v == NULL) \
        return NULL; \
\
    v->rows    = rows; \
    v->columns = 1; \
    v->min     = 1; \
    v->size    = rows; \
    v->type    = 0; \
    v->view    = 1; \
    v->M       = &A->M[i*rows]; \
\
    return v; \
}

/** @brief Get i-th column of matrix A
 *
 * @param [in] A matrix
 * @param [in] i column number
 * @retval x*y
 */
MATRIX_GET_COLUMN(matrix_get_column,         double,    matrix_t);

/** @brief Get i-th column of matrix A
 *
 * @param [in] A matrix
 * @param [in] i column number
 * @retval x*y
 */
MATRIX_GET_COLUMN(matrix_complex_get_column, complex_t, matrix_complex_t);


/** \defgroup sparse Functions for sparse matrices
 *  @{
 */

#ifdef SUPPORT_SPARSE

/** @brief Calculate eigenvalues of a sparse complex matrix
 *
 * @param [in] N     number of columns/rows of matrix
 * @param [in] nev   number of eigenvalues to compute
 * @param [in] which LM (largest magnitude), SM (smallest magnitude), LR (largest real part), SR (smallest real part), LI (largest imaginary part), SI (smallest imaginary part)
 * @param [in] Av callback function that implements the matrix-vector operation Av; the input vector is given as in, the vector Av must be written in out
 * @param [in,out] d on exit d contains the Rith approximations (must be of length nev+1)
 * @param [in] mxiter maximum number of Arnoldi update iterations allowed
 * @param [in] tol relative accuracy of the Ritz value
 * @param [in] data pointer that is given to callback function Av
 *
 * @retval 0 if successful
 */
int sparse_complex_eig(int N, int nev, char *which, void (*Av)(int N, complex_t *in, complex_t *out, void *data), complex_t *d, int mxiter, double tol, void *data)
{
    int ret = LIBHADES_ERROR_OOM;
    int info = 0;
    int ido = 0;
    char *bmat = "I";          /* standard eigenproblem */
    int ncv = MIN(2*nev+2, N);

    int ishift = 1;
    int mode = 1;
    int iparam[11] = { ishift, 0, mxiter, 1, 0, 0, mode, 0, 0, 0, 0 };

    int ipntr[14];
    int lworkl = ncv*(3*ncv + 5);
    int rvec = 0;
    char *howmny = "P";

    complex_t *workd = NULL, *workl = NULL, *resid = NULL, *v = NULL, *workev = NULL;
    int *select = NULL;
    double *rwork = NULL;

    /* allocate memory */
    ret = LIBHADES_ERROR_OOM;

    workd = malloc_cb(3*N*sizeof(complex_t));
    if(workd == NULL)
        goto out;
    workl = malloc_cb(lworkl*sizeof(complex_t));
    if(workl == NULL)
        goto out;
    rwork = malloc_cb(ncv*sizeof(double));
    if(rwork == NULL)
        goto out;
    resid  = malloc_cb(N*sizeof(complex_t));
    if(resid == NULL)
        goto out;
    v = malloc_cb(N*ncv*sizeof(complex_t));
    if(v == NULL)
        goto out;
    select = malloc_cb(ncv*sizeof(int));
    if(select == NULL)
        goto out;
    workev = malloc_cb((2*ncv)*sizeof(complex_t));
    if(workev == NULL)
        goto out;

    /* loop */
    while(1)
    {
        /* http://www.caam.rice.edu/software/ARPACK/UG/node138.html */
        znaupd_(
            &ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &N, iparam,
            ipntr, workd, workl, &lworkl, rwork, &info, strlen(bmat),
            strlen(which)
        );

        if(ido == 1 || ido == -1)
            Av(N, &workd[ipntr[0]-1], &workd[ipntr[1]-1], data);
        else if(ido == 99)
            break;
        else
        {
            ret = ido;
            goto out;
        }
    }

    if(info != 0)
    {
        ret = info;
        goto out;
    }

    /* http://www.mathkeisan.com/usersguide/man/zneupd.html */
    zneupd_(
        &rvec, howmny, select, d, NULL, &N, NULL, workev, bmat, &N, which,
        &nev, &tol, resid, &ncv, v, &N, iparam, ipntr, workd, workl, &lworkl,
        rwork, &info, strlen(howmny), strlen(bmat), strlen(which)
    );

    ret = info;

out:
    if(resid != NULL)
        free_cb(resid);
    if(v != NULL)
        free_cb(v);
    if(workd != NULL)
        free_cb(workd);
    if(workl != NULL)
        free_cb(workl);
    if(rwork != NULL)
        free_cb(rwork);
    if(select != NULL)
        free_cb(select);
    if(workev != NULL)
        free_cb(workev);

    return ret;
}

#endif
/** @}*/

/** \defgroup io Save/load matrices
 *  @{
 */

#define MATRIX_LOAD_FROM_STREAM(FUNCTION_NAME, MATRIX_TYPE, TYPE, ALLOC, TRANSPOSE, IS_COMPLEX) \
MATRIX_TYPE *FUNCTION_NAME(FILE *stream, int *ret) \
{ \
    MATRIX_TYPE *M; \
    uint16_t len; \
    int rows, columns, fortran_order, is_complex; \
    char header[10] = { 0 }; \
    char dict[2048] = { 0 }; \
\
    if(ret != NULL) \
        *ret = 0; \
\
    /* read magic string, major and minor number */ \
    fread(header, 8, 1, stream); \
    if(memcmp(header, "\x93NUMPY\x01\x00", 8) != 0) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_HEADER; \
        return NULL; \
    } \
\
    /* read length of dict */ \
    fread(&len, sizeof(uint16_t), 1, stream); \
\
    if(len >= sizeof(dict)/sizeof(dict[0])) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_INV_LENGTH; \
        return NULL; \
    } \
\
    fread(dict, sizeof(char), len, stream); \
\
    if(npy_dict_get_fortran_order(dict, &fortran_order) != 0) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_ORDER; \
        return NULL; \
    } \
\
    if(npy_dict_get_shape(dict, &rows, &columns) != 0) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_SHAPE; \
        return NULL; \
    } \
\
    if(npy_dict_get_descr(dict, &is_complex) != 0) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_DESCR; \
        return NULL; \
    } \
\
    if(is_complex != IS_COMPLEX) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_FORMAT; \
        return NULL; \
    } \
\
    M = ALLOC(rows,columns); \
    fread(M->M, sizeof(TYPE), rows*columns, stream); \
\
    if(!fortran_order) \
        TRANSPOSE(M); \
\
    M->min  = MIN(rows,columns); \
    M->size = rows*columns; \
    M->view = 0; \
    M->type = 0; \
\
    return M; \
}

/** @brief Load real matrix from stream
 *
 * Load real matrix A from a stream. This function will also allocate memory
 * for the matrix.
 *
 * If error != NULL, error will be set to:
 *  0                         if successful
 *  LIBHADES_ERROR_HEADER     if magic or major/minor number is invalid
 *  LIBHADES_ERROR_INV_LENGTH if length of dictionary is invalid (too long)
 *  LIBHADES_ERROR_ORDER      if order is invalid (Fortran/C order)
 *  LIBHADES_ERROR_SHAPE      if shape is invalid (rows/columns)
 *  LIBHADES_ERROR_DESCR      if dtype is wrong/not supported
 *  LIBHADES_ERROR_FORMAT     if wrong format (real instead of complex)
 *  0 rows/columns wrong
 *
 * @param [in] stream file handle of a opened file
 * @param [out] error error code
 * @retval A matrix
 * @retval NULL if an error occured
 */
MATRIX_LOAD_FROM_STREAM(matrix_load_from_stream, matrix_t, double, matrix_alloc, matrix_transpose, 0);

/** @brief Load complex matrix from stream
 *
 * Load complex matrix A from a stream. This function will also allocate memory
 * for the matrix.
 *
 * If error != NULL, error will be set to:
 *  0                         if successful
 *  LIBHADES_ERROR_HEADER     if magic or major/minor number is invalid
 *  LIBHADES_ERROR_INV_LENGTH if length of dictionary is invalid (too long)
 *  LIBHADES_ERROR_ORDER      if order is invalid (Fortran/C order)
 *  LIBHADES_ERROR_SHAPE      if shape is invalid (rows/columns)
 *  LIBHADES_ERROR_DESCR      if dtype is wrong/not supported
 *  LIBHADES_ERROR_FORMAT     if wrong format (complex instead of real)
 *  0 rows/columns wrong
 *
 * @param [in] stream file handle of a opened file
 * @param [out] error error code
 * @retval A matrix
 * @retval NULL if an error occured
 */
MATRIX_LOAD_FROM_STREAM(matrix_complex_load_from_stream, matrix_complex_t, complex_t, matrix_complex_alloc, matrix_complex_transpose, 1);


#define MATRIX_LOAD(FUNCTION_NAME, MATRIX_TYPE, LOAD_FUNCTION) \
MATRIX_TYPE *FUNCTION_NAME(const char *filename, int *ret) \
{ \
    FILE *stream; \
    MATRIX_TYPE *M; \
\
    if((stream = fopen(filename, "r")) == NULL) \
    { \
        if(ret != NULL) \
            *ret = LIBHADES_ERROR_IO; \
        return NULL; \
    } \
\
    M = LOAD_FUNCTION(stream, ret); \
\
    fclose(stream); \
\
    return M; \
}

/** @brief Load real matrix from file
 *
 * Load real matrix A from file given by filename. This function will also
 * allocate memory for the matrix. See \ref matrix_load_from_stream for errors.
 *
 * @param [in] filename path to the file
 * @param [out] error error code
 * @retval A matrix
 * @retval NULL if an error occured
 */
MATRIX_LOAD(matrix_load, matrix_t, matrix_load_from_stream);

/** @brief Load complex matrix from file
 *
 * Load complex matrix A from file given by filename. This function will also
 * allocate memory for the matrix. See \ref matrix_complex_load_from_stream for
 * errors.
 *
 * @param [in] filename path to the file
 * @param [out] error error code
 * @retval A matrix
 * @retval NULL if an error occured
 */
MATRIX_LOAD(matrix_complex_load, matrix_complex_t, matrix_complex_load_from_stream);

#define MATRIX_SAVE_TO_STREAM(FUNCTION_NAME, TYPE, MATRIX_TYPE, DTYPE) \
void FUNCTION_NAME(MATRIX_TYPE *M, FILE *stream) \
{ \
    char d_str[512] = { 0 }; \
    uint16_t len = 0; \
    const int rows = M->rows, columns = M->columns; \
\
    /* write magic string, major number and minor number */ \
    fwrite("\x93NUMPY\x01\x00", sizeof(char), 8, stream); \
\
    /* write length of header and header */ \
    snprintf(d_str, sizeof(d_str)/sizeof(d_str[0]), "{'descr': '%s', 'fortran_order': True, 'shape': (%d, %d), }", DTYPE, rows, columns); \
\
    len = strlen(d_str); \
\
    fwrite(&len,  sizeof(len),  1,   stream); \
    fwrite(d_str, sizeof(char), len, stream); \
\
    /* write matrix */ \
    fwrite(M->M, sizeof(TYPE), M->size, stream); \
}

/** @brief Save real matrix A to stream
 *
 * Save real matrix A to a file handle given by stream. The datatype
 * corresponds to Numpy's npy file format.
 *
 * @param [in] A      matrix to be dumped to file
 * @param [in] stream file handle of opened file
 */
MATRIX_SAVE_TO_STREAM(matrix_save_to_stream, double, matrix_t, "<d8");

/** @brief Save complex matrix A to stream
 *
 * Save complex matrix A to a file handle given by stream. The datatype
 * corresponds to Numpy's npy file format.
 *
 * @param [in] A      matrix to be dumped to file
 * @param [in] stream file handle of opened file
 */
MATRIX_SAVE_TO_STREAM(matrix_complex_save_to_stream, complex_t, matrix_complex_t, "<c16");


#define MATRIX_SAVE(FUNCTION_NAME, MATRIX_TYPE, SAVE_FUNCTION) \
int FUNCTION_NAME(MATRIX_TYPE *M, const char *filename) \
{ \
    FILE *stream = fopen(filename, "w"); \
    if(stream == NULL) \
        return LIBHADES_ERROR_IO; \
\
    SAVE_FUNCTION(M, stream); \
\
    fclose(stream); \
\
    return 0; \
};

/** @brief Save real matrix A to file
 *
 * Save real matrix A to file given by filename. The datatype corresponds to
 * Numpy's .npy file format.
 *
 * @param [in] A        matrix to be dumped to file
 * @param [in] filename path to the file
 * @retval 0 if successful
 * @retval LIBHADES_ERROR_IO if file could not be opened
 */
MATRIX_SAVE(matrix_save, matrix_t, matrix_save_to_stream);

/** @brief Save complex matrix A to file
 *
 * Save complex matrix A to file given by filename. The datatype corresponds to
 * Numpy's .npy file format.
 *
 * @param [in] A        matrix to be dumped to file
 * @param [in] filename path to the file
 * @retval 0 if successful
 * @retval LIBHADES_ERROR_IO if file could not be opened
 */
MATRIX_SAVE(matrix_complex_save, matrix_complex_t, matrix_complex_save_to_stream);

/** @}*/
