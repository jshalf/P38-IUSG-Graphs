#include "Jacobi.hpp"
#include "Matrix.hpp"

double JacobiRelax(CSR A, double *b, double **x, double *x_prev, int i);
double JacobiRelaxVolatile(CSR A, double *b, volatile double **x, volatile double *x_prev, int i);
double BlockJacobiRelax(CSR A, double *b, double **x, int i);
double BlockJacobiRelaxAtomic(CSR A, double *b, double **x, int i);
double BlockJacobiRelaxVolatile(CSR A, double *b, volatile double **x, int i);

void Jacobi(SolverData *solver, CSR A, double *b, double **x)
{
   int solver_type = solver->input.solver_type;
   int num_iters = solver->input.num_iters;
   int n = A.n;
   double *x_prev = (double *)calloc(n, sizeof(double));
   volatile double *x_vol = (double *)calloc(n, sizeof(double));
   volatile double *x_prev_vol = (double *)calloc(n, sizeof(double));

   #pragma omp parallel for
   for (int i = 0; i < n; i++){
      x_prev[i] = (*x)[i];
   }
   double start = omp_get_wtime();
   #pragma omp parallel
   {
      if (solver_type == SYNC_JACOBI){
         for (int iter = 0; iter < num_iters; iter++){
            #pragma omp for
            for (int i = 0; i < n; i++){
               (*x)[i] = JacobiRelax(A, b, x, x_prev, i);
            }
            #pragma omp for
            for (int i = 0; i < n; i++){
               x_prev[i] = (*x)[i];
            }
         }
      }
      if (solver_type == SYNC_BLOCK_JACOBI){

      }
      else if (solver_type == ASYNC_JACOBI){
         if (solver->input.atomic_flag){
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = JacobiRelaxVolatile(A, b, &x_vol, x_prev_vol, i);
                  #pragma omp atomic write
                  x_vol[i] = xi;
               }
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi;
                  xi = x_vol[i];
                  #pragma omp atomic write
                  x_prev_vol[i] = xi;
               }
            }
         }
         else {
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = JacobiRelax(A, b, x, x_prev, i);
                  (*x)[i] = xi;
               }
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  x_prev[i] = (*x)[i];
               }
            }
         }
      }
      else if (solver_type == ASYNC_BLOCK_JACOBI){
         if (solver->input.atomic_flag){
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = BlockJacobiRelaxVolatile(A, b, &x_vol, i);
                  #pragma omp atomic write
                  x_vol[i] = xi;
               }
            }
         }
         else {
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = BlockJacobiRelax(A, b, x, i);
                  (*x)[i] = xi;
               }
            }
         }
      }
   }
   solver->output.solve_wtime = omp_get_wtime() - start; 

   if (solver->input.atomic_flag &&
       (solver_type == ASYNC_JACOBI || solver_type == ASYNC_BLOCK_JACOBI)){
      for (int i = 0; i < n; i++){
         (*x)[i] = x_vol[i];
      }
   }

   free(x_prev);
   free((void *)x_vol);
   free((void *)x_prev_vol);
}

double JacobiRelax(CSR A, double *b, double **x, double *x_prev, int i)
{
   double res = b[i];
   for (int jj = A.i[i]; jj < A.i[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * x_prev[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double BlockJacobiRelax(CSR A, double *b, double **x, int i)
{
   double res = b[i];
   for (int jj = A.i[i]; jj < A.i[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * (*x)[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double JacobiRelaxVolatile(CSR A, double *b, volatile double **x, volatile double *x_prev, int i)
{
   double res = b[i];
   for (int jj = A.i[i]; jj < A.i[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * x_prev[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double BlockJacobiRelaxVolatile(CSR A, double *b, volatile double **x, int i)
{
   double res = b[i];
   for (int jj = A.i[i]; jj < A.i[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * (*x)[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double BlockJacobiRelaxAtomic(CSR A, double *b, volatile double **x, int i)
{
   double res = b[i];
   for (int jj = A.i[i]; jj < A.i[i+1]; jj++){
      int ii = A.j[jj];
      double xii;
      #pragma omp atomic read
      xii = (*x)[ii];
      res -= A.data[jj] * xii;
   }
   return (*x)[i] + res / A.diag[i];
}
