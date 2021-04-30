#include "Jacobi.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/***************************************************************
 * Solve for x in Ax=b using synchronous or asynchronous Jacobi.
 * A is in compressed sparse row (CSR) format.
 ***************************************************************/

double JacobiRelax(CSR A, double *b, double **x, double *x_prev, int i);
double AsyncJacobiRelax(CSR A, double *b, double **x, int i);
double AsyncJacobiRelaxAtomic(CSR A, double *b, double **x, int i);
double AsyncJacobiRelaxMsgQ(CSR A, double *b, double *x_ghost, Queue *Q, int i);

void Jacobi(SolverData *solver,
            CSR A, /* sparse matrix */
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

   int q_size;
   vector<vector<int>> put_targets(n);
   double *x_ghost;
   Queue Q;
   /* set up message queues */
   if (solver->input.MsgQ_flag == 1){
      x_ghost = (double *)malloc(nnz * sizeof(double));
      for (int i = 0; i < n; i++){
         for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
            int ii = A.j[jj];
            put_targets[ii].push_back(jj); /* put targets correspond to non-zeros in this row */
            x_ghost[jj] = (*x)[ii];
         }
      }
      q_size = nnz;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
   }

   double start = omp_get_wtime();
   #pragma omp parallel
   {
      if (solver_type == SYNC_JACOBI){ /* synchronous Jacobi */
         /* iterate until num_iters (naive convergence detection) */
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
      else if (solver_type == ASYNC_JACOBI){
         if (solver->input.MsgQ_flag == 1){ /* asynchronous Jacobi with message queues */
            /* iterate until num_iters (naive convergence detection) */
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = (*x)[i];
                  xi += AsyncJacobiRelaxMsgQ(A, b, x_ghost, &Q, i); /* relaxation of element i of x */
                  /* send element i of x using put primitive */
                  for (int j = 0; j < put_targets[i].size(); j++){
                     qPut(&Q, put_targets[i][j], xi);
                  }
                  (*x)[i] = xi;
               }
            }
         }
         else {
            if (solver->input.atomic_flag){ /* asynchronous Jacobi with atomics */
               /* iterate until num_iters (naive convergence detection) */
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double xi = AsyncJacobiRelaxAtomic(A, b, x, i); /* relaxation of element i of x */
                     /* atomically write to memory */
                     #pragma omp atomic write
                     (*x)[i] = xi;
                  }
               }
            }
            else { /* atomic asynchronous Jacobi without atomics */
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double xi = AsyncJacobiRelax(A, b, x, i); /* relaxation of element i of x */
                     (*x)[i] = xi; /* write to memeory without atomics */
                  }
               }
            }
         }
      }
   }
   solver->output.solve_wtime = omp_get_wtime() - start; 

   if (solver->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }

   if (solver_type == SYNC_JACOBI){
      free(x_prev);
   }
}

/* **********************************************************
 * relaxation routines for different Jacobi implementations
 * (See JacobiRelax() for detailed comments) 
 * **********************************************************/ 
double JacobiRelax(CSR A, /* sparse matrix */
                   double *b, /* right-hadn side */
                   double **x, /* current approximation to the solution (output) */
                   double *x_prev, /* approximation from previous iteration */
                   int i /* row to relax */
                   )
{
   /* compute residual for row i */
   double res = b[i]; /* initialize residual */
   for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){ /* loop over non-zeros in this row */
      int ii = A.j[jj]; /* column index */
      res -= A.data[jj] * x_prev[ii]; /* decrement residual */
   }
   /* return relaxation */
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelax(CSR A, double *b, double **x, int i)
{
   double res = b[i];
   for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * (*x)[ii];
   }
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelaxAtomic(CSR A, double *b, double **x, int i)
{
   double res = b[i];
   for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
      int ii = A.j[jj];
      double xii;
      /* atomically read element A.j[jj] of x */
      #pragma omp atomic read
      xii = (*x)[ii];
      res -= A.data[jj] * xii;
   }
   return (*x)[i] + res / A.diag[i];
}

double AsyncJacobiRelaxMsgQ(CSR A, double *b, double *x_ghost, Queue *Q, int i)
{
   double res = b[i];
   for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
      int ii = A.j[jj];
      double xii;
      qGet(Q, jj, &(x_ghost[jj])); /* get element A.j[jj] of x and save into x_ghost */
      res -= A.data[jj] * x_ghost[jj];
   }
   return res / A.diag[i];
}
