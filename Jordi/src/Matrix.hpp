#ifndef MATRIX_H
#define MATRIX_H

#define MATRIX_NONSYMMETRIC 0
#define MATRIX_LOWER 1
#define MATRIX_UPPER 2

#define PROBLEM_RANDOM 0
#define PROBLEM_5PT_POISSON 1

#include "Main.hpp"

typedef struct{
   int *i_ptr;
   int *i;
   int *j;
   double *data;
   double *diag;
   int n;
   int m;
   int nnz;
}CSR;

void Laplace_2D_5pt(InputData input, CSR *A, int n);

void RandomMatrix(InputData input, CSR *A, int n, int max_row_nnz, int mat_type);

void PrintCOO(CSR A, char *filename);

#endif
