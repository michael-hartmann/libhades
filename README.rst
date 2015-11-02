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


Examples
--------

Allocate a real 3x3 matrix:
.. code:: C
matrix_t *M = matrix_alloc(3,3);


Create a complex 4x4 matrix and initialize it as identity matrix:
.. code:: C
matrix_t *M = matrix_complex_eye(3,3);


Multiply two complex matrices and create a new matrix C=A*B:
.. code:: C
matrix_complex_t *C = matrix_complex_mult(A,B,NULL)


Calculate eigenvalues and right eigenvectors of a generic complex matrix M:
.. code:: C
eig_complex_generic(M, w, NULL, vr);


Installation
------------


Documentation
-------------


How to contribute
-----------------



Bug reports
-----------


License information
-------------------
