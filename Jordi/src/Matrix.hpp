#ifndef MATRIX_H
#define MATRIX_H

/* matrix pattern types */
#define MATRIX_NONSYMMETRIC 0
#define MATRIX_LOWER 1
#define MATRIX_UPPER 2

/* test matrix types */
#define PROBLEM_RANDOM 0
#define PROBLEM_5PT_POISSON 1
#define PROBLEM_FILE 2

#include "Main.hpp"

/* CSR struct */
typedef struct{
   int *start; /* pointer to row starts */
   int *i; /* columns indices (used in COO and CSC) */
   int *j; /* row indices (used in COO and CSR) */
   double *data; /* matrix values */
   double *diag; /* diagonal elements */
   int n; /* number of rows */
   int m; /* number of columns */
   int nnz; /* number of non-zero values */
}Matrix;

void Laplace_2D_5pt(InputData input, Matrix *A, int n);

void RandomMatrix(InputData input,
                  Matrix *A,
                  int n,
                  int max_row_nnz,
                  int mat_type,
                  int csc_flag,
                  int coo_flag);

void PrintMatrix(Matrix A,
                 char *filename,
                 int print_diag_flag,
                 int csc_flag);

void freadBinaryMatrix(char *mat_file_str,
                       Matrix *A,
                       int include_diag_flag,
                       int csc_flag,
                       int coo_flag,
                       int mat_type);

#endif
