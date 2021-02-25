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
   int *i_ptr; /* pointer to row starts (only used in CSR algorithms) */
   int *i; /* row indices (only used in COO algorithms) */
   int *j; /* columns indices */
   double *data; /* matrix values */
   double *diag; /* diagonal elements */
   int n; /* number of rows */
   int m; /* number of columns */
   int nnz; /* number of non-zero values */
}CSR;

void Laplace_2D_5pt(InputData input, CSR *A, int n);

void RandomMatrix(InputData input, CSR *A, int n, int max_row_nnz, int mat_type);

void PrintCOO(CSR A, char *filename, int print_diag_flag);

void freadBinaryMatrix(char *mat_file_str, CSR *A, int include_diag_flag);

#endif
