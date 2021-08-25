#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <vector>

// Fortran function declarations
extern "C" {
    extern int dgetrf_(int*, int*, double*, int*, int*, int*);
    extern int dgetri_(int*, double*, int*, int*, double*, int*, int*);
    extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
}

// C BLAS/LAPACK wrappers
int C_DGETRF(int m, int n, double* a, int lda, int* ipiv) {
    int info;
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int C_DGETRI(int n, double* a, int lda, int* ipiv, double* work, int lwork) {
    int info;
    dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha, double* a, int lda, double* b,
             int ldb, double beta, double* c, int ldc) {
    if (m == 0 || n == 0 || k == 0) return;
    dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
}

// Simply matrix multiply using BLAS
std::vector<double> matrix_multiply(std::vector<double> L, std::vector<double> R, int dim) {
    assert(L.size() == R.size());
    std::vector<double> product(L);
    C_DGEMM('n', 'n', dim, dim, dim, 1.0, L.data(), dim, R.data(), dim, 0.0, product.data(), dim);
    return product;
}

// A copy of Psi4's general_invert() routine
std::vector<double> invert(const std::vector<double> &matrix, int dim) {
    int lwork = dim * dim;
    std::vector<double> work(lwork);
    std::vector<int> ipiv(dim);
    std::vector<double> inverse(matrix);

    int err = C_DGETRF(dim, dim, inverse.data(), dim, ipiv.data());
    if (err != 0) {
        if (err < 0) {
            printf("invert: C_DGETRF: argument %d has invalid parameter.\n", -err);
            abort();
        }
        if (err > 1) {
            printf(
                "invert: C_DGETRF: the (%d,%d) element of the factor U or L is "
                "zero, and the inverse could not be computed.\n",
                err, err);
            abort();
        }
    }

    err = C_DGETRI(dim, inverse.data(), dim, ipiv.data(), work.data(), lwork);
    if (err != 0) {
        if (err < 0) {
            printf("invert: C_DGETRI: argument %d has invalid parameter.\n", -err);
            abort();
        }
        if (err > 1) {
            printf(
                "invert: C_DGETRI: the (%d,%d) element of the factor U or L is "
                "zero, and the inverse could not be computed.\n",
                err, err);
            abort();
        }
    }
    return inverse;
}

void print_matrix(const char *header, int dim, const std::vector<double> &matrix, int nrows, int ncols) {
    printf("%s\n", header);
    for (int row = 0; row < nrows; ++row) {
        for (int col = 0; col < ncols; ++col) {
            printf("%20.10f ", matrix[row * dim + col]);
        }
        printf("\n");
    }
}

int main() {
    std::ifstream matrix_file("matrix.bin", std::ios::in | std::ios::binary);
    if (matrix_file) {
        // Get file size
        matrix_file.seekg(0, std::ios::end);
        size_t file_size = matrix_file.tellg();
        int N2 = file_size / sizeof(double);
        int N = std::sqrt(N2);
        matrix_file.seekg(0, std::ios::beg);
        printf("Matrix is %d x %d\n", N, N);

        // Read the matrix
        matrix_file.seekg(0, std::ios::end);
        std::vector<double> matrix(N2);
        matrix_file.seekg(0, std::ios::beg);
        matrix_file.read(reinterpret_cast<char*>(matrix.data()), file_size);
        matrix_file.close();

        print_matrix("Input matrix", N, matrix, 8, 8);
        auto inverse = invert(matrix, N);
        print_matrix("Inverted matrix", N, inverse, 8, 8);
        auto identity = matrix_multiply(matrix, inverse, N);
        print_matrix("Matrix x Inverted matrix", N, identity, 8, 8);

        double min_diagonal = 999;
        double max_diagonal = -999;
        double min_offdiagonal = 999;
        double max_offdiagonal = -999;
        for (int row = 0; row < N; ++row) {
            min_diagonal = std::min(min_diagonal, identity[row * N + row]);
            max_diagonal = std::max(max_diagonal, identity[row * N + row]);
            for (int col = row+1; col < N; ++col) {
                min_offdiagonal = std::min(min_offdiagonal, identity[row * N + col]);
                max_offdiagonal = std::max(max_offdiagonal, identity[row * N + col]);
            }
        }
        printf("Diagonal elements:\n");
        printf("Max: %20.10f Min: %20.10f\n", max_diagonal, min_diagonal);
        printf("Off-diagonal elements:\n");
        printf("Max: %20.10f Min: %20.10f\n", max_offdiagonal, min_offdiagonal);
    } else {
        printf("matrix.bin does not exist, or is unreadable\n");
        abort();
    }
    return 0;
}
