#include "Jacobi.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/***************************************************************
 * Solve for x in Ax=b using synchronous or asynchronous Jacobi.
 * A is in compressed sparse row (Matrix) format.
 ***************************************************************/

double JacobiRelax_CSR(Matrix A, double *b, double **x, double *x_prev, int i);
double AsyncJacobiRelax_CSR(Matrix A, double *b, double **x, int i);
double AsyncJacobiRelaxAtomic_CSR(Matrix A, double *b, double **x, int i);
double AsyncJacobiRelaxMsgQ_CSR(Matrix A, double *b, double *x_ghost, Queue *Q, int i);

void JacobiRelaxAtomic_CSC(Matrix A, double **r, double z, int i);
void AsyncJacobiRelax_CSC(Matrix A, double **r, double z, int i);

void Jacobi(SolverData *solver,
            Matrix A, /* sparse matrix */
            double *b, /* right-hand side */
            double **x /* solution (output) */
            )
{
   int solver_type = solver->input.solver_type;
   int num_iters = solver->input.num_iters;
   int n = A.n;
   int nnz = A.nnz;

   /* for synchronous Jacobi, x from previous iteration must be saved */
   double *x_prev;
   if (solver_type == SYNC_JACOBI){
      x_prev = (double *)calloc(n, sizeof(double));
      #pragma omp parallel for
      for (int i = 0; i < n; i++){
         x_prev[i] = (*x)[i];
      }
   }

   double *r;
   if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
      r = (double *)calloc(n, sizeof(double));
   }

   int q_size;
   vector<vector<int>> put_targets(n);
   double *x_ghost;
   Queue Q;
   /* set up message queues */
   if (solver->input.MsgQ_flag == 1){
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSR){
         x_ghost = (double *)malloc(nnz * sizeof(double));
         for (int i = 0; i < n; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               int ii = A.j[jj];
               put_targets[ii].push_back(jj); /* put targesolver correspond to non-zeros in this row */
               x_ghost[jj] = (*x)[ii];
            }
         }
      }
      q_size = nnz;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         #pragma omp for
         for (int i = 0; i < n; i++){
            r[i] = b[i];
         }
      }

      double solve_start = omp_get_wtime();   
      if (solver_type == SYNC_JACOBI){ /* synchronous Jacobi */
         if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for
               for (int i = 0; i < n; i++){
                  double z;
                  #pragma omp atomic read
                  z = r[i];
                  z /= A.diag[i];
                  (*x)[i] += z;
                  JacobiRelaxAtomic_CSC(A, &r, z, i);
               }
            }
         }
         else {
            /* iterate until num_iters (naive convergence detection) */
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for
               for (int i = 0; i < n; i++){
                  (*x)[i] = JacobiRelax_CSR(A, b, x, x_prev, i);
               }
               #pragma omp for
               for (int i = 0; i < n; i++){
                  x_prev[i] = (*x)[i];
               }
            }
         }
      }
      else if (solver_type == ASYNC_JACOBI){
         if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
            if (solver->input.MsgQ_flag == 1){ /* asynchronous Jacobi with message queues */
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double z;
                     while (qGet(&Q, i, &z)){
                        r[i] -= z;
                     }
                     z = r[i];
                     z /= A.diag[i];
                     (*x)[i] += z;
                     for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                        qPut(&Q, A.i[jj], A.data[jj] * z);
                     }
                  }
               }      
            }
            else {
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double z;
                     if (solver->input.atomic_flag == 1){
                        #pragma omp atomic read
                        z = r[i];
                     }
                     else {
                        z = r[i];
                     }
                     z /= A.diag[i];
                     (*x)[i] += z;
                     if (solver->input.atomic_flag == 1){
                        JacobiRelaxAtomic_CSC(A, &r, z, i);
                     }
                     else {
                        AsyncJacobiRelax_CSC(A, &r, z, i);
                     }
                  }
               }
            }
         }
         else {
            if (solver->input.MsgQ_flag == 1){ /* asynchronous Jacobi with message queues */
               /* iterate until num_iters (naive convergence detection) */
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double xi = (*x)[i];
                     xi += AsyncJacobiRelaxMsgQ_CSR(A, b, x_ghost, &Q, i); /* relaxation of element i of x */
                     /* send element i of x using put primitive */
                     for (int j = 0; j < put_targets[i].size(); j++){
                        qPut(&Q, put_targets[i][j], xi);
                     }
                     (*x)[i] = xi;
                  }
               }
            }
            else {
               /* iterate until num_iters (naive convergence detection) */
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double xi;
                     if (solver->input.atomic_flag){ /* atomically write to memory */
                        xi = AsyncJacobiRelaxAtomic_CSR(A, b, x, i); /* relaxation of element i of x */
                        #pragma omp atomic write
                        (*x)[i] = xi;
                     }
                     else {
                        xi = AsyncJacobiRelax_CSR(A, b, x, i); /* relaxation of element i of x */
                        (*x)[i] = xi;
                     }
                  }
               }
            }
         }
      }
      solver->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;
   }

   if (solver->input.MsgQ_flag == 1){
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSR){
         free(x_ghost);
      }
      qDestroyLock(&Q);
      qFree(&Q);
   }

   if (solver_type == SYNC_JACOBI){
      free(x_prev);
   }

   if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
      free(r);
   }
}

/* **********************************************************
 * relaxation routines for different Jacobi implementations
 * (See JacobiRelax() for detailed commensolver) 
 * **********************************************************/ 
double JacobiRelax_CSR(Matrix A, /* sparse matrix */
                       double *b, /* right-hadn side */
                       double **x, /* current approximation to the solution (output) */
                       double *x_prev, /* approximation from previous iteration */
                       int i /* row to relax */
                       )
{
   /* compute residual for row i */
   double res = b[i]; /* initialize residual */
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){ /* loop over non-zeros in this row */
      int ii = A.j[jj]; /* column index */
      res -= A.data[jj] * x_prev[ii]; /* decrement residual */
   }
   /* return relaxation */
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelaxAtomic_CSR(Matrix A, double *b, double **x, int i)
{
   double res = b[i];
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.j[jj];
      double xii;
      /* atomically read element A.j[jj] of x */
      #pragma omp atomic read
      xii = (*x)[ii];
      res -= A.data[jj] * xii;
   }
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelax_CSR(Matrix A, double *b, double **x, int i)
{
   double res = b[i];
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * (*x)[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelaxMsgQ_CSR(Matrix A, double *b, double *x_ghost, Queue *Q, int i)
{
   double res = b[i];
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.j[jj];
      double xii;
      qGet(Q, jj, &(x_ghost[jj])); /* get element A.j[jj] of x and save into x_ghost */
      res -= A.data[jj] * x_ghost[jj];
   }
   return res / A.diag[i];
}

void JacobiRelaxAtomic_CSC(Matrix A, double **r, double z, int i)
{
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.i[jj];
      #pragma omp atomic
      (*r)[ii] -= A.data[jj] * z;
   }
}

void AsyncJacobiRelax_CSC(Matrix A, double **r, double z, int i)
{
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.i[jj];
      (*r)[ii] -= A.data[jj] * z;
   }
}
