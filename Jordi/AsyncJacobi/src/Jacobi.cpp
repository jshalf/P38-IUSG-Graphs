#include "Jacobi.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

/***************************************************************
 * Solve for x in Ax=b using synchronous or asynchronous Jacobi.
 * A is in compressed sparse row (Matrix) format.
 ***************************************************************/

double JacobiRelax_CSR(SparseMatrix A, double *b, double **x, double *x_prev, int i);
double AsyncJacobiRelax_CSR(SparseMatrix A, double *b, double **x, int i);
double AsyncJacobiRelaxAtomic_CSR(SparseMatrix A, double *b, double **x, int i);

double AsyncJacobiRelaxMsgQ_CSR(SparseMatrix A, double *b, double *x_ghost, MessageQueue<double> *Q, int i);

void JacobiRelaxAtomic_CSC(SparseMatrix A, double **r, double z, int i);
void JacobiRelax_CSC(SparseMatrix A, double **r, double z, int i);

void Jacobi(SolverData *solver,
            SparseMatrix A, /* sparse matrix */
            double *b, /* right-hand side */
            double **x /* solution (output) */
            )
{
   int solver_type = solver->input.solver_type;
   int num_iters = solver->input.num_iters;
   int lump = 1;

   int n = A.GetNumRows();
   int nnz = A.GetNNZ();

   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<int> row_idx = A.GetRowIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();

   /* for synchronous Jacobi, x from previous iteration must be saved */
   double *x_prev;
   //if (solver_type == SYNC_JACOBI){
      x_prev = (double *)calloc(n, sizeof(double));
      #pragma omp parallel for
      for (int i = 0; i < n; i++){
         x_prev[i] = (*x)[i];
      }
   //}

   double *r = (double *)calloc(n, sizeof(double));
   if (A.GetStorageType() == SparseMatrixStorageType::CSC){
      #pragma omp parallel for
      for (int i = 0; i < n; i++){
         r[i] = b[i];
      }
   }

   double *x_ghost;
   MessageQueue<double> *Q, *Q_comm;
   /* set up message queues */
   if (solver->input.MsgQ_flag == 1){
      if (A.GetStorageType() == SparseMatrixStorageType::CSC){
         Q = new MessageQueue<double>(n);
      }
      else {
         if (A.GetMatType() == MatrixType::symm){
            Q = new MessageQueue<double>(n);
         }
         else {
            x_ghost = (double *)malloc(nnz * sizeof(double));
            for (int i = 0; i < n; i++){
               for (int jj = start[i]; jj < start[i+1]; jj++){
                  int ii = col_idx[jj];
                  x_ghost[jj] = (*x)[ii];
               }
            }
            Q = new MessageQueue<double>(nnz);
            Q_comm = new MessageQueue<double>(n);
         }
      }
   }

   #pragma omp parallel
   {
      double comp_wtime_start, comp_wtime_stop, comp_wtime = 0.0;
      double MsgQ_wtime_start, MsgQ_wtime_stop, MsgQ_wtime = 0.0;
      double MsgQ_put_wtime_start, MsgQ_put_wtime_stop, MsgQ_put_wtime = 0.0;
      double MsgQ_get_wtime_start, MsgQ_get_wtime_stop, MsgQ_get_wtime = 0.0;
      double solve_start, solve_stop;
      uint64_t MsgQ_cycles_start, MsgQ_cycles_stop, MsgQ_cycles = 0;
      uint64_t MsgQ_put_cycles_start, MsgQ_put_cycles_stop, MsgQ_put_cycles = 0;
      uint64_t MsgQ_get_cycles_start, MsgQ_get_cycles_stop, MsgQ_get_cycles = 0;
      double dummy = 0.0;
      int num_qPuts = 0, num_qGets = 0;
      int dummy_num_qPuts = 0, dummy_num_qGets = 0, dummy_num_spins = 0, dummy_mean_num_spins;
      int i_loc, jj_loc, n_loc, nnz_loc;
      double *r_loc, *x_loc, *z_loc;
      int *my_rows;

      int tid = omp_get_thread_num();

      n_loc = 0, nnz_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         n_loc++;
         for (int jj = start[i]; jj < start[i+1]; jj++){
            nnz_loc++;
         }
      }

      my_rows = (int *)calloc(n_loc, sizeof(int));
      r_loc = (double *)calloc(n_loc, sizeof(double));

      i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         my_rows[i_loc] = i;
         r_loc[i_loc] = b[i];
         i_loc++;
      }

      double *x_prev_loc = (double *)calloc(nnz_loc, sizeof(double));

      //A_loc.diag = (double *)calloc(n_loc, sizeof(double));
      //A_loc.start = (int *)calloc(n_loc+1, sizeof(int));
      //A_loc.data = (double *)calloc(nnz_loc, sizeof(double));
      //if (A.GetStorageType() == SparseMatrixStorageType::CSC){
      //   A_loc.i = (int *)calloc(nnz_loc, sizeof(int));
      //}
      //else {
      //   A_loc.j = (int *)calloc(nnz_loc, sizeof(int));
      //} 
      //i_loc = 0, jj_loc = 0;
      //#pragma omp for schedule(static, lump) nowait
      //for (int i = 0; i < n; i++){
      //   A_loc.diag[i_loc] = diag[i];
      //   int i_nnz = 0;
      //   for (int jj = start[i]; jj < start[i+1]; jj++){
      //      if (A.GetStorageType() == SparseMatrixStorageType::CSC){
      //         A_loc.i[jj_loc] = row_idx[jj];
      //      }
      //      else {
      //         A_loc.j[jj_loc] = col_idx[jj];
      //      }
      //      A_loc.data[jj_loc] = mat_values[jj];
      //      jj_loc++;
      //      i_nnz++;
      //   }
      //   A_loc.start[i_loc+1] = A_loc.start[i_loc] + i_nnz;
      //   i_loc++;
      //}

      if (solver->input.MsgQ_flag == 1){
         z_loc = (double *)calloc(nnz_loc, sizeof(double));
         x_loc = (double *)calloc(n_loc, sizeof(double));
      }

      solve_start = omp_get_wtime();   


/*************************
 * 
 *       SYNC JACOBI
 *
 *************************/
      if (solver_type == SYNC_JACOBI){
         if (A.GetStorageType() == SparseMatrixStorageType::CSC){
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for
               for (int i = 0; i < n; i++){
                  double z;
                  #pragma omp atomic read
                  z = r[i];
                  z /= diag[i];
                  (*x)[i] += z;
                  JacobiRelaxAtomic_CSC(A, &r, z, i);
               }
            }
         }
         else {
            /* iterate until num_iters (naive convergence detection) */
            for (int iter = 0; iter < num_iters; iter++){
               #pragma omp for schedule(static, lump)
               for (int i = 0; i < n; i++){
                  (*x)[i] = JacobiRelax_CSR(A, b, x, x_prev, i);
               }
               #pragma omp for schedule(static, lump)
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
         if (A.GetStorageType() == SparseMatrixStorageType::CSC){
            /*************************
             * 
             * Async Jacobi MsgQ CSC
             *
             *************************/
            if (solver->input.MsgQ_flag == 1){
               for (int iter = 0; iter < num_iters; iter++){
                  #pragma omp for schedule(static, lump) nowait
                  for (int i = 0; i < n; i++){
                     double z;
                     int get_flag = Q->qGet(i, &z);
                     while (get_flag == 1){
                        r[i] -= z;
                        get_flag = Q->qGet(i, &z);
                     }
                     z = r[i];
                     z /= diag[i];
                     (*x)[i] += z;
                     for (int jj = start[i]; jj < start[i+1]; jj++){
                        Q->qPut(row_idx[jj], mat_values[jj] * z);
                     }
                  }
               }
            }










            /***************************
             * 
             * Async Jacobi atomic CSC
             *
             ***************************/
            else {
               if (solver->input.atomic_flag == 1){ /* atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        #pragma omp atomic read
                        z = r[i];
                        z /= diag[i];
                        (*x)[i] += z;
                        JacobiRelaxAtomic_CSC(A, &r, z, i);
                     }
                  }
               }
               else { /* no atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double z;
                        z = r[i];
                        z /= diag[i];
                        (*x)[i] += z;
                        JacobiRelax_CSC(A, &r, z, i);
                     }
                  }
               }
            }
         }










         else {
            /***********************************
             * 
             * SYMMETRIC Async Jacobi MsgQ CSR
             *
             ***********************************/
            if (solver->input.MsgQ_flag == 1){
               if (A.GetMatType() == MatrixType::symm){
                  /*****************
                   *   MsgQ wtime
                   *****************/
                  if (solver->input.MsgQ_wtime_flag == 1){
                     MsgQ_wtime_start = omp_get_wtime();
                     int mean_num_spins = 0;
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z = r_loc[i_loc] / diag[i];
                           x_loc[i_loc] += z;
                           for (int jj = start[i]; jj < start[i+1]; jj++){
                              int j = col_idx[jj];
                              double y = mat_values[jj] * z;
                              if (i == j){
                                 r_loc[i_loc] -= y;
                              }
                              else {
                                 Q->qPut(j, y);
                                 num_qPuts++;
                              }
                           }

                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           int num_spins = 0;
                           while (1){
                              get_flag = Q->qGet(i, &z_recv);
                              if (!get_flag) break;
                              z_accum += z_recv;
                              num_qGets++;
                              num_spins++;
                           }
                           mean_num_spins += num_spins;
                           r_loc[i_loc] -= z_accum;
                        }
                     }
                     mean_num_spins /= n_loc * num_iters;
                     MsgQ_wtime_stop = omp_get_wtime();
                     MsgQ_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;

                     MsgQ_wtime_start = omp_get_wtime();
                     dummy_mean_num_spins = 0;
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z = r_loc[i_loc] / diag[i];
                           x_loc[i_loc] += z;
                           for (int jj = start[i]; jj < start[i+1]; jj++){
                              int j = col_idx[jj];
                              double y = mat_values[jj] * z;
                              if (i == j){
                                 r_loc[i_loc] -= y;
                              }
                              else {
                                 dummy += j + y;
                                 dummy_num_qPuts++;
                              }
                           }

                           int get_flag;
                           double z_accum = 0.0, z_recv = 1.0;
                           int dummy_num_spins = 0;
                           for (int s = 0; s < mean_num_spins; s++){
                              z_accum += z_recv;
                              dummy_num_qGets++;
                              dummy_num_spins++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }
                     }
                     dummy_mean_num_spins /= n_loc * num_iters;
                     MsgQ_wtime_stop = omp_get_wtime();
                     MsgQ_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;
                  }
                  /*****************
                   *   MsgQ cycles
                   *****************/
                  else if (solver->input.MsgQ_cycles_flag == 1){
                  }
                  /*****************
                   *   comp wtime
                   *****************/
                  else if (solver->input.comp_wtime_flag == 1){
                     int mean_num_spins = 0;
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z = r_loc[i_loc] / diag[i];
                           x_loc[i_loc] += z;
                           for (int jj = start[i]; jj < start[i+1]; jj++){
                              int j = col_idx[jj];
                              double y = mat_values[jj] * z;
                              if (i == j){
                                 r_loc[i_loc] -= y;
                              }
                              else {
                                 Q->qPut(j, y);
                                 num_qPuts++;
                              }
                           }

                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           int num_spins = 0;
                           while (1){
                              get_flag = Q->qGet(i, &z_recv);
                              if (!get_flag) break;
                              z_accum += z_recv;
                              num_qGets++;
                              num_spins++;
                           }
                           mean_num_spins += num_spins;
                           r_loc[i_loc] -= z_accum;
                        }
                     }
                     mean_num_spins /= n_loc * num_iters;

                     comp_wtime_start = omp_get_wtime();
                     dummy_mean_num_spins = 0;
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z = r_loc[i_loc] / diag[i];
                           x_loc[i_loc] += z;
                           for (int jj = start[i]; jj < start[i+1]; jj++){
                              int j = col_idx[jj];
                              double y = mat_values[jj] * z;
                              if (i == j){
                                 r_loc[i_loc] -= y;
                              }
                              else {
                                 dummy += j + y;
                                 dummy_num_qPuts++;
                              }
                           }

                           int get_flag;
                           double z_accum = 0.0, z_recv = 1.0;
                           int dummy_num_spins = 0;
                           for (int s = 0; s < mean_num_spins; s++){
                              z_accum += z_recv;
                              dummy_num_qGets++;
                              dummy_num_spins++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }
                     }
                     dummy_mean_num_spins /= n_loc * num_iters;
                     comp_wtime_stop = omp_get_wtime();
                     comp_wtime += comp_wtime_stop - comp_wtime_start;
                  }
                  /******************************
                   *   MsgQ no-op and comp no-op
                   ******************************/
                  else if ((solver->input.MsgQ_noop_flag == 1) && (solver->input.comp_noop_flag == 1)){
                  }
                  /*****************
                   *   MsgQ no-op
                   *****************/
                  else if (solver->input.MsgQ_noop_flag == 1){
                  }
                  /*****************
                   *   comp no-op
                   *****************/
                  else if (solver->input.comp_noop_flag == 1){
                  }
                  /**********************************************
                   * standard scheme (no timers, no-ops, etc...)
                   **********************************************/
                  else {
                     int mean_num_spins = 0;
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z = r_loc[i_loc] / diag[i];
                           x_loc[i_loc] += z; 
                           for (int jj = start[i]; jj < start[i+1]; jj++){
                              int j = col_idx[jj];
                              double y = mat_values[jj] * z;
                              if (i == j){
                                 r_loc[i_loc] -= y;
                              }
                              else {
                                 Q->qPut(j, y);
                                 num_qPuts++;
                              }
                           }
                           
                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           int num_spins = 0;
                           while (1){
                              get_flag = Q->qGet(i, &z_recv);
                              if (!get_flag) break;
                              z_accum += z_recv;
                              num_qGets++;
                              num_spins++;
                           }
                           mean_num_spins += num_spins;
                           r_loc[i_loc] -= z_accum;
                        }
                     }
                     mean_num_spins /= n_loc * num_iters;
                  }
               }










               /**************************************
                *
                * NON-SYMMETRIC Async Jacobi MsgQ CSR
                *
                * ************************************/
               else {
                  int n_loc = 0;
                  #pragma omp for schedule(static, lump) nowait
                  for (int i = 0; i < n; i++) n_loc++;
                  vector<vector<int>> put_targets(n_loc);
                  #pragma omp for schedule(static, lump) nowait
                  for (int i = 0; i < n; i++){
                     for (int jj = start[i]; jj < start[i+1]; jj++){
                        int ii = col_idx[jj];
                        Q_comm->qPut(ii, (double)jj);
                     }
                  }
                  /* iterate until num_iters (naive convergence detection) */
                  for (int iter = 0; iter < num_iters; iter++){
                     int i_loc = 0;
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double target;
                        if (Q_comm->qGet(i, &target)){
                           put_targets[i_loc].push_back((int)target);
                        }
                        double xi = (*x)[i];
                        xi += AsyncJacobiRelaxMsgQ_CSR(A, b, x_ghost, Q, i); /* relaxation of element i of x */
                        /* send element i of x using put primitive */
                        for (int j = 0; j < put_targets[i_loc].size(); j++){
                           Q->qPut(put_targets[i_loc][j], xi);
                        }
                        (*x)[i] = xi;
                        i_loc++;
                     }
                  }
               }
            }











            /**************************
             * 
             * Async Jacobi atomic CSR
             *
             **************************/
            else {
               if (solver->input.atomic_flag == 1){ /* atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double res = b[i];
                        for (int jj = start[i]; jj < start[i+1]; jj++){
                           int ii = col_idx[jj];
                           double xii;
                           /* atomically read element col_idx[jj] of x */
                           #pragma omp atomic read
                           xii = (*x)[ii];
                           res -= mat_values[jj] * xii;
                        }
                        #pragma omp atomic
                        (*x)[i] += res / diag[i];
                     }
                  }
               }
               else { /* no atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double res = b[i]; /* initialize residual */
                        for (int jj = start[i]; jj < start[i+1]; jj++){ /* loop over non-zeros in this row */
                           int ii = col_idx[jj]; /* column index */
                           res -= mat_values[jj] * (*x)[ii]; /* decrement residual */
                        }
                        (*x)[i] += res / diag[i];
                     }
                  }
               }
            }
         }
      }
      solve_stop = omp_get_wtime();

      solver->output.solve_wtime_vec[tid] = solve_stop - solve_start;

      if (solver->input.MsgQ_flag == 1){
         for (i_loc = 0; i_loc < n_loc; i_loc++){
            int i = my_rows[i_loc];
            (*x)[i] = x_loc[i_loc];
         }
         MsgQ_put_wtime = MsgQ_wtime;

         //if (solver->input.MsgQ_wtime_flag == 1){
            solver->output.MsgQ_put_wtime_vec[tid] = MsgQ_put_wtime;
            solver->output.MsgQ_get_wtime_vec[tid] = MsgQ_get_wtime;
         //}
         //else if (solver->input.MsgQ_cycles_flag == 1){
            solver->output.MsgQ_put_cycles_vec[tid] = MsgQ_put_cycles;
            solver->output.MsgQ_get_cycles_vec[tid] = MsgQ_get_cycles;
         //}
         //else if (solver->input.comp_wtime_flag == 1){
            solver->output.comp_wtime_vec[tid] = comp_wtime;
         //}

         solver->output.num_qGets_vec[tid] = num_qGets;
         solver->output.num_qPuts_vec[tid] = num_qPuts;

         free(z_loc);
         free(x_loc);

         dummy += (dummy_num_qGets + dummy_num_qPuts + dummy_num_spins + dummy_mean_num_spins);

         PrintDummy(dummy);

         #pragma omp barrier
      }

      #pragma omp barrier
      jj_loc = 0;
      #pragma omp for schedule(static, lump) nowait
       for (int i = 0; i < n; i++){
          for (int jj = start[i]; jj < start[i+1]; jj++){
             int ii = col_idx[jj];
             (*x)[ii] += x_prev_loc[jj_loc];
             jj_loc++;
          }
       }

      free(x_prev_loc);
      free(my_rows);
      free(r_loc);
      //free(A_loc.diag);
      //free(A_loc.data);
      //free(A_loc.start);
      //if (A.GetStorageType() == SparseMatrixStorageType::CSC){
      //   free(A_loc.i);
      //}
      //else {
      //   free(A_loc.j);
      //}
   }

   if (solver->input.MsgQ_flag == 1){
      //if (A.GetStorageType() == MATRIX_STORAGE_CSR){
      //   if (A.GetMatType() == MatrixType::symm){
      //   }
      //   else {
      //      free(x_ghost);
      //      qDestroyLock(&Q_comm);
      //      qFree(&Q_comm);
      //   }
      //}
      //qDestroyLock(&Q);
      //qFree(&Q);
   }

   //if (solver_type == SYNC_JACOBI){
      free(x_prev);
   //}

   free(r);
}

/*****************************************************************************
 * relaxation routines for different Jacobi implementations
 * (the comments in JacobiRelax_CSR() can be extended to all functions below) 
 *****************************************************************************/ 
double JacobiRelax_CSR(SparseMatrix A, /* sparse matrix */
                       double *b, /* right-hadn side */
                       double **x, /* current approximation to the solution (output) */
                       double *x_prev, /* approximation from previous iteration */
                       int i /* row to relax */
                       )
{
   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();

   /* compute residual for row i */
   double res = b[i]; /* initialize residual */
   for (int jj = start[i]; jj < start[i+1]; jj++){ /* loop over non-zeros in this row */
      int ii = col_idx[jj]; /* column index */
      res -= mat_values[jj] * x_prev[ii]; /* decrement residual */
   }
   /* return relaxation */
   return (*x)[i] + res / diag[i];
}

double AsyncJacobiRelaxAtomic_CSR(SparseMatrix A, double *b, double **x, int i)
{
   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();

   double res = b[i];
   for (int jj = start[i]; jj < start[i+1]; jj++){
      int ii = col_idx[jj];
      double xii;
      /* atomically read element col_idx[jj] of x */
      //#pragma omp atomic read
      xii = (*x)[ii];
      res -= mat_values[jj] * xii;
   }
   return (*x)[i] + res / diag[i];
}

double AsyncJacobiRelax_CSR(SparseMatrix A, double *b, double **x, int i)
{
   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();

   double res = b[i];
   for (int jj = start[i]; jj < start[i+1]; jj++){
      int ii = col_idx[jj];
      res -= mat_values[jj] * (*x)[ii];
   }
   return (*x)[i] + res / diag[i];
}

double AsyncJacobiRelaxMsgQ_CSR(SparseMatrix A, double *b, double *x_ghost, MessageQueue<double> *Q, int i)
{
   vector<int> start = A.GetIndexStarts();
   vector<int> col_idx = A.GetColIndices();
   vector<double> mat_values = A.GetValues();
   vector<double> diag = A.GetDiagValues();

   double res = b[i];
   for (int jj = start[i]; jj < start[i+1]; jj++){
      int ii = col_idx[jj];
      double xii;
      Q->qGet(jj, &(x_ghost[jj])); /* get element col_idx[jj] of x and save into x_ghost */
      res -= mat_values[jj] * x_ghost[jj];
   }
   return res / diag[i];
}

void JacobiRelaxAtomic_CSC(SparseMatrix A, double **r, double z, int i)
{
   vector<int> start = A.GetIndexStarts();
   vector<int> row_idx = A.GetRowIndices();
   vector<double> mat_values = A.GetValues();

   for (int jj = start[i]; jj < start[i+1]; jj++){
      int ii = row_idx[jj];
      #pragma omp atomic
      (*r)[ii] -= mat_values[jj] * z;
   }
}

void JacobiRelax_CSC(SparseMatrix A, double **r, double z, int i)
{
   vector<int> start = A.GetIndexStarts();
   vector<int> row_idx = A.GetRowIndices();
   vector<double> mat_values = A.GetValues();

   for (int jj = start[i]; jj < start[i+1]; jj++){
      int ii = row_idx[jj];
      (*r)[ii] -= mat_values[jj] * z;
   }
}
