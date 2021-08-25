About
=====

This is a minimal test case to attempt to recreate a [specific
problem](https://github.com/psi4/psi4/pull/2283) diagonalizing a matrix in a
Psi4 test case.

Compiling
=========

Intel Compilers
---------------
```bash
icpc -std=c++14 -mkl test.cc
```

GCC with LAPACK
---------------
```bash
g++ test.cc  -std=c++14  /your/path/to/lapack/lib/liblapack.a /your/path/to/lapack/lib/libblas.a /your/path/to/gcc/lib64/libgfortran.so
```
(the last part is only needed if BLAS/LAPACK need libgfortran)

macOS
-----
```bash
g++ test.cc  -std=c++14 -framework Accelerate
```
