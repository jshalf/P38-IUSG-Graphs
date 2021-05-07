#include "Jacobi.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

/***************************************************************
 * Solve for x in Ax=b using synchronous or asynchronous Jacobi.
 * A is in compressed sparse row (Matrix) format.
 ***************************************************************/

double JacobiRelax_CSR(Matrix A, double *b, double **x, double *x_prev, int i);
double AsyncJacobiRelax_CSR(Matrix A, double *b, double **x, int i);
double AsyncJacobiRelaxAtomic_CSR(Matrix A, double *b, double **x, int i);
double AsyncJacobiRelaxMsgQ_CSR(Matrix A, double *b, double *x_ghost, Queue *Q, int i);

void JacobiRelaxAtomic_CSC(Matrix A, double **r, double z, int i);
void JacobiRelax_CSC(Matrix A, double **r, double z, int i);

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
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         q_size = n;
      }
      else {
         x_ghost = (double *)malloc(nnz * sizeof(double));
         for (int i = 0; i < n; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               int ii = A.j[jj];
               put_targets[ii].push_back(jj); /* put targets correspond to non-zeros in this row */
               x_ghost[jj] = (*x)[ii];
            }
         }
         q_size = nnz;
      }
      qAlloc(&Q, q_size);
      qInitLock(&Q);
   }

   #pragma omp parallel
   {
      double comp_wtime_start, comp_wtime = 0.0, MsgQ_wtime_start, MsgQ_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles = 0, comp_cycles_start, comp_cycles = 0;

      int tid = omp_get_thread_num();
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         #pragma omp for
         for (int i = 0; i < n; i++){
            r[i] = b[i];
         }
      }

      double solve_start = omp_get_wtime();   


      /*************************
       * 
       *       SYNC JACOBI
       *
       *************************/
      if (solver_type == SYNC_JACOBI){
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






/*****************************************************
 *
 *                       ASYNC JACOBI
 *
 *****************************************************/
      else if (solver_type == ASYNC_JACOBI){
         if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){



            /*************************
             * Async Jacobi MsgQ CSC
             *************************/
            if (solver->input.MsgQ_flag == 1){
               /*****************
                *   MsgQ wtime
                *****************/
               if (solver->input.MsgQ_wtime_flag == 1){
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        MsgQ_wtime_start = omp_get_wtime();
                        int get_flag = qGet(&Q, i, &z);
                        MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                        while (get_flag == 1){
                           r[i] -= z;
                           MsgQ_wtime_start = omp_get_wtime();
                           get_flag = qGet(&Q, i, &z);
                           MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                        }
                        z = r[i];
                        z /= A.diag[i];
                        (*x)[i] += z;
                        for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                           MsgQ_wtime_start = omp_get_wtime();
                           qPut(&Q, A.i[jj], A.data[jj] * z);
                           MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                        }
                     }
                  }
               }
               /*****************
                *   MsgQ cycles
                *****************/
               else if (solver->input.MsgQ_cycles_flag == 1){
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        MsgQ_cycles_start = rdtsc();
                        int get_flag = qGet(&Q, i, &z);
                        MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                        while (get_flag == 1){
                           r[i] -= z;
                           MsgQ_cycles_start = rdtsc();
                           get_flag = qGet(&Q, i, &z);
                           MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                        }
                        z = r[i];
                        z /= A.diag[i];
                        (*x)[i] += z;
                        for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                           MsgQ_cycles_start = rdtsc();
                           qPut(&Q, A.i[jj], A.data[jj] * z);
                           MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                        }
                     }
                  }
               }
               /*****************
                *   comp wtime
                *****************/
               else if (solver->input.comp_wtime_flag == 1){
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        int get_flag = qGet(&Q, i, &z);
                        while (get_flag == 1){
                           comp_wtime_start = omp_get_wtime();
                           r[i] -= z;
                           comp_wtime += omp_get_wtime() - comp_wtime_start;
                           get_flag = qGet(&Q, i, &z);
                        }
                        comp_wtime_start = omp_get_wtime();
                        z = r[i];
                        z /= A.diag[i];
                        (*x)[i] += z;
                        comp_wtime += omp_get_wtime() - comp_wtime_start;
                        for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                           qPut(&Q, A.i[jj], A.data[jj] * z);
                        }
                     }
                  }
               }
               /*****************
                *   MsgQ no-op
                *****************/
               else if (solver->input.MsgQ_noop_flag == 1){
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z = 0.0;
                        r[i] -= z;
                        z = r[i];
                        z /= A.diag[i];
                        (*x)[i] += z;
                        for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                        }
                     }
                  }
               }
               /*****************
                *   comp no-op
                *****************/
               else if (solver->input.comp_noop_flag == 1){
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        int get_flag = qGet(&Q, i, &z);
                        while (get_flag == 1){
                           get_flag = qGet(&Q, i, &z);
                        }
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
                        int get_flag = qGet(&Q, i, &z);
                        while (get_flag == 1){
                           r[i] -= z;
                           get_flag = qGet(&Q, i, &z);
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
            }

            /***************************
             * Async Jacobi atomic CSC
             ***************************/
            else {
               if (solver->input.atomic_flag == 1){ /* atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
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
               else { /* no atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        z = r[i];
                        z /= A.diag[i];
                        (*x)[i] += z;
                        JacobiRelax_CSC(A, &r, z, i);
                     }
                  }
               }
            }
         }


         else {
            /*************************
             * Async Jacobi MsgQ CSR
             *************************/
            if (solver->input.MsgQ_flag == 1){
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
            /**************************
             * Async Jacobi atomic CSR
             **************************/
            else {
               if (solver->input.atomic_flag == 1){ /* atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double xi;
                        xi = AsyncJacobiRelaxAtomic_CSR(A, b, x, i); /* relaxation of element i of x */
                        #pragma omp atomic write
                        (*x)[i] = xi;
                     }
                  }
               }
               else { /* no atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for nowait
                     for (int i = 0; i < n; i++){
                        double xi;
                        xi = AsyncJacobiRelax_CSR(A, b, x, i); /* relaxation of element i of x */
                        (*x)[i] = xi;
                     }
                  }
               }
            }
         }
      }

      solver->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      if (solver->input.MsgQ_flag == 1){
         if (solver->input.MsgQ_wtime_flag == 1){
            solver->output.MsgQ_wtime_vec[tid] = MsgQ_wtime;
         }
         else if (solver->input.MsgQ_cycles_flag == 1){
            solver->output.MsgQ_cycles_vec[tid] = MsgQ_cycles;
         }
         else if (solver->input.comp_wtime_flag == 1){
            solver->output.comp_wtime_vec[tid] = comp_wtime;
         }
         else if (solver->input.comp_cycles_flag == 1){
            solver->output.comp_cycles_vec[tid] = comp_cycles;
         }
      }
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

void JacobiRelax_CSC(Matrix A, double **r, double z, int i)
{
   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
      int ii = A.i[jj];
      (*r)[ii] -= A.data[jj] * z;
   }
}
