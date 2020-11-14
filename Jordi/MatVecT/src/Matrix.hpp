#ifndef MATRIX_H
#define MATRIX_H

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

void Laplace_2D_5pt(MatVecData *mv, CSR *A, int n);

#endif
