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

double JacobiRelax_CSR_loc(Matrix A_loc, double *b_loc, double **x_loc, double *x_prev, int i, int i_loc);
double AsyncJacobiRelax_CSR_loc(Matrix A_loc, double *b_loc, double **x, int i, int i_loc);
double AsyncJacobiRelaxAtomic_CSR_loc(Matrix A_loc, double *b_loc, double **x, int i, int i_loc);

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
   int lump = 1;

   /* for synchronous Jacobi, x from previous iteration must be saved */
   double *x_prev;
   if (solver_type == SYNC_JACOBI){
      x_prev = (double *)calloc(n, sizeof(double));
      #pragma omp parallel for
      for (int i = 0; i < n; i++){
         x_prev[i] = (*x)[i];
      }
   }

   double *r = (double *)calloc(n, sizeof(double));
   if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
      #pragma omp parallel for
      for (int i = 0; i < n; i++){
         r[i] = b[i];
      }
   }

   double *x_ghost;
   Queue Q, Q_comm;
   /* set up message queues */
   if (solver->input.MsgQ_flag == 1){
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         qAlloc(&Q, n);
         qInitLock(&Q);
      }
      else {
         if (solver->input.symm_flag == 1){
            qAlloc(&Q, n);
            qInitLock(&Q);
         }
         else {
            x_ghost = (double *)malloc(nnz * sizeof(double));
            for (int i = 0; i < n; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  int ii = A.j[jj];
                  x_ghost[jj] = (*x)[ii];
               }
            }
            qAlloc(&Q_comm, n);
            qInitLock(&Q_comm);
            qAlloc(&Q, nnz);
            qInitLock(&Q);
         }
      }
   }

   #pragma omp parallel
   {
      double comp_wtime_start, comp_wtime_stop, comp_wtime = 0.0;
      double MsgQ_wtime_start, MsgQ_wtime_stop, MsgQ_wtime = 0.0;
      double MsgQ_put_wtime_start, MsgQ_put_wtime_stop, MsgQ_put_wtime = 0.0;
      double MsgQ_get_wtime_start, MsgQ_get_wtime_stop, MsgQ_get_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles_stop, MsgQ_cycles = 0;
      uint64_t MsgQ_put_cycles_start, MsgQ_put_cycles_stop, MsgQ_put_cycles = 0;
      uint64_t MsgQ_get_cycles_start, MsgQ_get_cycles_stop, MsgQ_get_cycles = 0;
      double dummy = 0.0;
      int num_qPuts = 0, num_qGets = 0;
      int dummy_num_qPuts = 0, dummy_num_qGets = 0;
      int i_loc, jj_loc, n_loc, nnz_loc;
      double *r_loc, *x_loc, *z_loc;
      int *my_rows;
      Matrix A_loc;

      int tid = omp_get_thread_num();

      n_loc = 0, nnz_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         n_loc++;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            nnz_loc++;
         }
      }

      my_rows = (int *)calloc(n_loc, sizeof(int));
      r_loc = (double *)calloc(n_loc, sizeof(double));
      A_loc.diag = (double *)calloc(n_loc, sizeof(double));
      A_loc.start = (int *)calloc(n_loc+1, sizeof(int));
      A_loc.data = (double *)calloc(nnz_loc, sizeof(double));

      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         A_loc.i = (int *)calloc(nnz_loc, sizeof(int));
      }
      else {
         A_loc.j = (int *)calloc(nnz_loc, sizeof(int));
      }
      
      i_loc = 0, jj_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         my_rows[i_loc] = i;
         r_loc[i_loc] = b[i];
         A_loc.diag[i_loc] = A.diag[i];
         int i_nnz = 0;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
               A_loc.i[jj_loc] = A.i[jj];
            }
            else {
               A_loc.j[jj_loc] = A.j[jj];
            }
            A_loc.data[jj_loc] = A.data[jj];
            jj_loc++;
            i_nnz++;
         }
         A_loc.start[i_loc+1] = A_loc.start[i_loc] + i_nnz;
         i_loc++;
      }

      if (solver->input.MsgQ_flag == 1){
         z_loc = (double *)calloc(nnz_loc, sizeof(double));
         x_loc = (double *)calloc(n_loc, sizeof(double));
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
         if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
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
                        z /= A.diag[i];
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
                        z /= A.diag[i];
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
               if (solver->input.symm_flag == 1){
                  /*****************
                   *   MsgQ wtime
                   *****************/
                  if (solver->input.MsgQ_wtime_flag == 1){
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }

                        MsgQ_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              qPut(&Q, j, z_loc[jj_loc]);
                              num_qPuts++;
                           }
                        }
                        MsgQ_wtime_stop = omp_get_wtime();
                        MsgQ_put_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;

                        MsgQ_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           get_flag = qGet(&Q, i, &z_recv);
                           if (get_flag == 1){
                              z_accum += z_recv;
                              num_qGets++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }
                        MsgQ_wtime_stop = omp_get_wtime();
                        MsgQ_get_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }

                        MsgQ_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              dummy += (j + z_loc[jj_loc]);
                              dummy_num_qPuts++;
                           }
                        }
                        MsgQ_wtime_stop = omp_get_wtime();
                        MsgQ_put_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;

                        MsgQ_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z_accum = 1.0, z_recv = 1.0;
                           z_accum += z_recv;
                           r_loc[i_loc] -= z_accum;
                           dummy_num_qGets++;
                           dummy += i;
                        }
                        MsgQ_wtime_stop = omp_get_wtime();
                        MsgQ_get_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;
                     }
                  }
                  /*****************
                   *   MsgQ cycles
                   *****************/
                  else if (solver->input.MsgQ_cycles_flag == 1){
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }

                        MsgQ_cycles_start = rdtsc();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              qPut(&Q, j, z_loc[jj_loc]);
                              num_qPuts++;
                           }
                        }
                        MsgQ_cycles_stop = rdtsc();
                        MsgQ_put_cycles += MsgQ_cycles_stop - MsgQ_cycles_start;

                        MsgQ_cycles_start = rdtsc();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           get_flag = qGet(&Q, i, &z_recv);
                           if (get_flag == 1){
                              z_accum += z_recv;
                              num_qGets++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }
                        MsgQ_cycles_stop = rdtsc();
                        MsgQ_get_cycles += MsgQ_cycles_stop - MsgQ_cycles_start;

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }

                        MsgQ_cycles_start = rdtsc();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              dummy += (j + z_loc[jj_loc]);
                              dummy_num_qPuts++;
                           }
                        }
                        MsgQ_cycles_stop = rdtsc(); 
                        MsgQ_put_cycles -= MsgQ_cycles_stop - MsgQ_cycles_start;

                        MsgQ_cycles_start = rdtsc();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z_accum = 1.0, z_recv = 1.0;
                           z_accum += z_recv;
                           r_loc[i_loc] -= z_accum;
                           dummy_num_qGets++;
                           dummy += i;
                        }
                        MsgQ_cycles_stop = rdtsc();
                        MsgQ_get_cycles -= MsgQ_cycles_stop - MsgQ_cycles_start;
                     }
                  }
                  /*****************
                   *   comp wtime
                   *****************/
                  else if (solver->input.comp_wtime_flag == 1){
                     for (int iter = 0; iter < num_iters; iter++){
                        comp_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }
                        comp_wtime_stop = omp_get_wtime();
                        comp_wtime += comp_wtime_stop - comp_wtime_start;

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              qPut(&Q, j, z_loc[jj_loc]);
                              num_qPuts++;
                           }
                        }

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           get_flag = qGet(&Q, i, &z_recv);
                           if (get_flag == 1){
                              z_accum += z_recv;
                              num_qGets++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z;
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }

                        comp_wtime_start = omp_get_wtime();
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              dummy += (j + z_loc[jj_loc]);
                              dummy_num_qPuts++;
                           }
                        }

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           double z_accum = 1.0;
                           r_loc[i_loc] -= z_accum;
                           dummy_num_qGets++;
                           dummy += i;
                        }
                        comp_wtime_stop = omp_get_wtime();
                        comp_wtime += comp_wtime_stop - comp_wtime_start;
                     }
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
                     for (int iter = 0; iter < num_iters; iter++){
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           double z = r_loc[i_loc] / A_loc.diag[i_loc];
                           x_loc[i_loc] += z; 
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              z_loc[jj_loc] = A_loc.data[jj_loc] * z;
                           }
                        }
                        
                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                              int j = A_loc.j[jj_loc];
                              qPut(&Q, j, z_loc[jj_loc]);
                              num_qPuts++;
                           }
                        }

                        for (i_loc = 0; i_loc < n_loc; i_loc++){
                           int i = my_rows[i_loc];
                           int get_flag;
                           double z_accum = 0.0, z_recv = 0.0;
                           get_flag = qGet(&Q, i, &z_recv);
                           if (get_flag == 1){
                              z_accum += z_recv;
                              num_qGets++;
                           }
                           r_loc[i_loc] -= z_accum;
                        }
                     }
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
                     for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                        int ii = A.j[jj];
                        qPut(&Q_comm, ii, (double)jj);
                     }
                  }
                  /* iterate until num_iters (naive convergence detection) */
                  for (int iter = 0; iter < num_iters; iter++){
                     int i_loc = 0;
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < n; i++){
                        double target;
                        if (qGet(&Q_comm, i, &target)){
                           put_targets[i_loc].push_back((int)target);
                        }
                        double xi = (*x)[i];
                        xi += AsyncJacobiRelaxMsgQ_CSR(A, b, x_ghost, &Q, i); /* relaxation of element i of x */
                        /* send element i of x using put primitive */
                        for (int j = 0; j < put_targets[i_loc].size(); j++){
                           qPut(&Q, put_targets[i_loc][j], xi);
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
                        double xi;
                        xi = AsyncJacobiRelaxAtomic_CSR(A, b, x, i); /* relaxation of element i of x */
                        #pragma omp atomic write
                        (*x)[i] = xi;
                     }
                  }
               }
               else { /* no atomic */
                  for (int iter = 0; iter < num_iters; iter++){
                     #pragma omp for schedule(static, lump) nowait
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
         for (i_loc = 0; i_loc < n_loc; i_loc++){
            int i = my_rows[i_loc];
            (*x)[i] = x_loc[i_loc];
         }

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

         dummy += (dummy_num_qGets + dummy_num_qPuts);

         PrintDummy(dummy);

         #pragma omp barrier
      }

      free(my_rows);
      free(r_loc);
      free(A_loc.diag);
      free(A_loc.data);
      free(A_loc.start);
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSC){
         free(A_loc.i);
      }
      else {
         free(A_loc.j);
      }
   }

   if (solver->input.MsgQ_flag == 1){
      if (solver->input.mat_storage_type == MATRIX_STORAGE_CSR){
         if (solver->input.symm_flag == 1){
         }
         else {
            free(x_ghost);
            qDestroyLock(&Q_comm);
            qFree(&Q_comm);
         }
      }
      qDestroyLock(&Q);
      qFree(&Q);
   }

   if (solver_type == SYNC_JACOBI){
      free(x_prev);
   }

   free(r);
}

/************************************************************
 * relaxation routines for different Jacobi implementations
 * (See JacobiRelax() for detailed commensolver) 
 ************************************************************/ 
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

double JacobiRelax_CSR_loc(Matrix A_loc, /* sparse matrix */
                           double *b_loc, /* right-hadn side */
                           double **x_loc, /* current approximation to the solution (output) */
                           double *x_prev, /* approximation from previous iteration */
                           int i, /* row to relax */
                           int i_loc
                           )
{
   /* compute residual for row i */
   double res = b_loc[i_loc]; /* initialize residual */
   for (int jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){ /* loop over non-zeros in this row */
      int ii = A_loc.j[jj_loc]; /* column index */
      res -= A_loc.data[jj_loc] * x_prev[ii]; /* decrement residual */
   }
   /* return relaxation */
   return (*x_loc)[i_loc] + res / A_loc.diag[i_loc];
}

double AsyncJacobiRelaxAtomic_CSR_loc(Matrix A_loc, double *b_loc, double **x, int i, int i_loc)
{
   double res = b_loc[i_loc];
   for (int jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
      int ii = A_loc.j[jj_loc];
      double xii;
      /* atomically read element A.j[jj] of x */
      #pragma omp atomic read
      xii = (*x)[ii];
      res -= A_loc.data[jj_loc] * xii;
   }
   return (*x)[i] + res / A_loc.diag[i_loc];
}

double AsyncJacobiRelax_CSR_loc(Matrix A_loc, double *b_loc, double **x, int i, int i_loc)
{
   double res = b_loc[i_loc];
   for (int jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
      int ii = A_loc.j[jj_loc];
      res -= A_loc.data[jj_loc] * (*x)[ii];
   }
   return (*x)[i] + res / A_loc.diag[i_loc];
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
