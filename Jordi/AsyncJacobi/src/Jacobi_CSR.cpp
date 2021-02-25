#include "Jacobi.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

double JacobiRelax(CSR A, double *b, double **x, double *x_prev, int i);
double AsyncJacobiRelax(CSR A, double *b, double **x, int i);
double AsyncJacobiRelaxAtomic(CSR A, double *b, double **x, int i);
double AsyncJacobiRelaxMsgQ(CSR A, double *b, double *x_ghost, Queue *Q, int i);

void Jacobi(SolverData *solver, CSR A, double *b, double **x)
{
   int solver_type = solver->input.solver_type;
   int num_iters = solver->input.num_iters;
   int n = A.n;
   int nnz = A.nnz;
   double *x_prev = (double *)calloc(n, sizeof(double));
   //volatile double *x_vol = (double *)calloc(n, sizeof(double));
   //volatile double *x_prev_vol = (double *)calloc(n, sizeof(double));

   #pragma omp parallel for
   for (int i = 0; i < n; i++){
      x_prev[i] = (*x)[i];
   }

   int q_size;
   vector<vector<int>> put_map(n);
   double *x_ghost;
   Queue Q;
   if (solver->input.MsgQ_flag == 1){
      x_ghost = (double *)malloc(nnz * sizeof(double));
      for (int i = 0; i < n; i++){
         for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
            int ii = A.j[jj];
            put_map[ii].push_back(jj);
         }
      }
      q_size = nnz;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
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
      else if (solver_type == ASYNC_JACOBI){
         if (solver->input.MsgQ_flag == 1){
            #pragma omp for
            for (int i = 0; i < n; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  int ii = A.j[jj];
                  x_ghost[jj] = (*x)[ii];
               }
            }
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for nowait
               for (int i = 0; i < n; i++){
                  double xi = (*x)[i];
                  xi += BlockJacobiRelaxMsgQ(A, b, x_ghost, &Q, i); 
                  for (int j = 0; j < put_map[i].size(); j++){
                     qPut(&Q, put_map[i][j], xi);
                  }
                  (*x)[i] = xi;
               }
            }
         }
         else {
            if (solver->input.atomic_flag){
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for nowait
                  for (int i = 0; i < n; i++){
                     double xi = AsyncJacobiRelaxAtomic(A, b, x, i);
                     #pragma omp atomic write
                     (*x)[i] = xi;
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
      if (solver->input.MsgQ_flag == 1){
         #pragma omp barrier
         #pragma omp for
         for (int i = 0; i < n; i++){
            double xi = 0;
            qGet(&Q, i, &xi);
            (*x)[i] = xi;
         }
      }
   }
   solver->output.solve_wtime = omp_get_wtime() - start; 

   if (solver->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }

   free(x_prev);
}

double JacobiRelax(CSR A, double *b, double **x, double *x_prev, int i)
{
   double res = b[i];
   for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
      int ii = A.j[jj];
      res -= A.data[jj] * x_prev[ii];
   }
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
      if(qPoll(Q, jj, &xii)){
         qGet(Q, jj, &(x_ghost[jj]));
      };
      res -= A.data[jj] * x_ghost[jj];
   }
   return res / A.diag[i];
}
