#ifndef MATRIX_H
#define MATRIX_H

#include "Main.hpp"

typedef struct{
   int *i;
   int *j;
   double *data;
   double *diag;
   int n;
   int nnz;
}CSR;

void Laplace_2D_5pt(CSR *A, int n);

#endif
