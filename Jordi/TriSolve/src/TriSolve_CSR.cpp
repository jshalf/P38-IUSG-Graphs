#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"

/* Serial TriSolve (sparse matrix must be in compressed sparse row (CSR) format) */
void TriSolve_CSR(TriSolveData *ts,
                  CSR L, /* lower triangular part of matrix */
                  CSR U, /* upper triangular part of matrix */
                  int *L_perm, /* ordering for computing elements of x */
                  int *U_perm, /* ordering for computing elements of y */
                  double *x, /* result of Lx=b (output) */
                  double *y, /* result of Uy=x (output) */
                  double *b /* right-hand side */
                  )
{
   int num_rows = L.n;

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

/* Level-scheduled TriSolve (sparse matrix must be in compressed sparse row (CSR) format) */
void TriSolve_LevelSets_CSR(TriSolveData *ts,
                            CSR L, /* lower triangular part of matrix */
                            CSR U, /* upper triangular part of matrix */
                            double *x, /* result of Lx=b (output) */
                            double *y, /* result of Uy=x (output) */
                            double *b)
{
   int num_rows = L.n;
   int lump = 1;

   #pragma omp parallel
   {
      for (int l = 0; l < ts->L_lvl_set.num_levels; l++){ /* loop over level sets L */
         #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
         for (int ii = ts->L_lvl_set.level_start[l]; ii < ts->L_lvl_set.level_start[l+1]; ii++){
            int i = ts->L_lvl_set.perm[ii];
            x[i] = b[i] / L.diag[i];
            for (int jj = L.i_ptr[i]; jj < L.i_ptr[i+1]; jj++){
               x[i] -= L.data[jj] * x[L.j[jj]] / L.diag[i];
            }
         }
      }
      for (int l = 0; l < ts->U_lvl_set.num_levels; l++){ /* loop over level sets U */
         #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
         for (int ii = ts->U_lvl_set.level_start[l]; ii < ts->U_lvl_set.level_start[l+1]; ii++){
            int i = ts->U_lvl_set.perm[ii];
            y[i] = x[i] / U.diag[i];
            for (int jj = U.i_ptr[i]; jj < U.i_ptr[i+1]; jj++){
               y[i] -= U.data[jj] * y[U.j[jj]] / U.diag[i];
            }
         }
      }
   }
}
