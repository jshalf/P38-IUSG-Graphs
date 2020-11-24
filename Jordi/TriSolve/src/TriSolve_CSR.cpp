#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"

void TriSolve_CSR(TriSolveData *ts,
                  CSR L,
                  CSR U,
                  int *L_perm,
                  int *U_perm,
                  double *x,
                  double *y,
                  double *b)
{
   int num_rows = L.n;

   for (int i = 0; i < num_rows; i++){
      x[i] = 0;
      y[i] = 0;
   }

   for (int ii = 0; ii < num_rows; ii++){
      int i = L_perm[ii];
      x[i] = b[i] / L.diag[i];
      for (int jj = L.i_ptr[i]; jj < L.i_ptr[i+1]; jj++){
         x[i] -= L.data[jj] * x[L.j[jj]] / L.diag[i];
      }
   }
   for (int ii = 0; ii < num_rows; ii++){
      int i = U_perm[ii];
      y[i] = x[i] / U.diag[i];
      for (int jj = U.i_ptr[i]; jj < U.i_ptr[i+1]; jj++){
         y[i] -= U.data[jj] * y[U.j[jj]] / U.diag[i];
      }
   }
}
