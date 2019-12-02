========
libhades
========

.. contents::

What is libhades?
-----------------

libhades is a small wrapper around LAPACK and BLAS. While LAPACK and
BLAS are very fast and reliable, the API is from the 70ies. For example,
zgeev, the routine to compute eigenvalues and left and/or right eigenvectors
of a generic complex matrix, has 14 arguments!

libhades tries to make things easier for you. It defines a simple matrix type
for complex and real matrices and gives you simple functions to add or multiply
matrices, to find eigenvalues and eigenvectors, calculate the matrix exponential
for arbitrary matrices, to solve differential equations, to invert matrices or
to calculate a LU decomposition and much more. While in the background still
LAPACK or BLAS do all the job, you have a simple and clean API. And in case you
really need to call a LAPACK function yourself, you can still do so.


Features
--------

- low-level operations on matrices like addition, multiplication...
- calculation of trace, determinant and norm of matrices
- LU decomposition and inverting of matrices
- Kronecker product
- eigenvalue problems
- matrix exponential using squaring and scaling algorithm with Padé-approximation
- integration of ordinary differential equations
- symplectic implicit integrator


Examples
--------

Allocate a real 3x3 matrix::

    matrix_t *M = matrix_alloc(3,3);

Create a complex 4x4 identity matrix::

    matrix_t *M = matrix_complex_eye(4,4);

Multiply two complex matrices and create a new matrix C=A*B::

    matrix_complex_t *C = matrix_complex_mult(A,B,NULL)

Calculate eigenvalues and right eigenvectors of a generic complex matrix M::

    eig_complex_generic(M, w, NULL, vr);

A full example to invert a real matrix (this is an example in the examples/
directory)::

    #include <stdio.h>
    #include <libhades.h>

    int main(int argc, char *argv[])
    {
        matrix_t *M = matrix_alloc(2,2);

        matrix_set(M, 0,0, 1);
        matrix_set(M, 0,1, 2);
        matrix_set(M, 1,0, 3);
        matrix_set(M, 1,1, 4);

        printf("Inverse of\n");
        matrix_fprintf(stdout, M, "%+4g", "  ", "\n");

        matrix_invert(M);
        
        printf("is:\n");
        matrix_fprintf(stdout, M, "%+4g", "  ", "\n");
         
        return 0;
    }

You can compile this file using::

    gcc -O2 -Wall invert.c -o invert -lm -lhades -llapack -lblas

There are more examples available in the examples/ directory. Just have a look!


Accessing LAPACK function
-------------------------

It's easy to access LAPACK functions. The matrix type is just a struct with
some information on the matrix, for example matrix_t is defined as::
    typedef struct {
        int rows;    /**< rows of matrix */
        int columns; /**< columns of matrix */
        int min;     /**< MIN(rows,columns) */
        int size;    /**< size = rows*columns */
        int type;    /**< type of matrix */
        int view;    /**< matrix is a view */
        double *M;   /**< data in column major order */
    } matrix_t;
 
M is a pointer the memory. The matrix is saved in column-major order (or
Fortran order). For this reason you don't have to transpose matrices. We store
the matrix just like LAPACK. To access matrix elements, use the macros
matrix_get and matrix set::

    M10 = matrix_get(M, 1,0); /* get matrix element in 2nd row, 1st column /*
    matrix_set(M, 2,3, 5);    /* set matrix element in 2nd row, 3rd column to 5 */

To call LAPACK function, just call the function name of the LAPACK function and
add an underscore at the end. So, if you want to call dgetrf, call dgetrf_.
There is one example available, so have a look at examples/lu_lapack.c.


Installation
------------

You need the development version of LAPACK and BLAS installed on your computer. On
Ubuntu/Debian you can install the dependencies using::

    $ apt-get install gcc libc6-dev make libblas-dev liblapack-dev

At the moment there is no build system. To compile the library change to the
directory libhades/ and run::

    make
    make install

This will compile the library and copy the shared object file libhades.so to
/usr/lib. You might need to run ldconfig afterwards.


Documentation
-------------

Documentation is available using Doxygen at
https://michael-hartmann.github.io/libhades/html/.


How to contribute
-----------------

Send bug reports, feature requests and merge requests! libhades is still in
development and there are probably a lot of bugs. I'm also happy for unit
tests.


License information
-------------------

The code is in the public domain, see the LICENSE file.


Bibliography
------------

- Moler, Loan, "Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later", SIAM Review, 2005
- Awad H. Al-Mohy and Nicholas J. Higham (2009) "A New Scaling and Squaring Algorithm for the Matrix Exponential." SIAM Journal on Matrix Analysis and Applications. 31 (3). pp. 970-989. ISSN 1095-7162
- Higham, "Functions of Matrices: Theory and Computation", Society for Industrial and Applied Mathematics, 2008
- Markiewicz, "Survey On Symplectic Integrators"
