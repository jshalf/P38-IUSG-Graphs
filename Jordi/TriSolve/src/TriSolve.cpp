#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/******************
 * Serial TriSolve 
 ******************/
void TriSolve_Seq(TriSolveData *ts,
                  Matrix T, /* triangular matrix */
                  int *T_perm, /* ordering for computing elements of x */
                  double *x, /* solution (output) */
                  double *b /* right-hand side */
                  )
{
   int n = T.n;

   if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
      for (int i = 0; i < n; i++) x[i] = b[i];
      for (int I = 0; I < n; I++){
         int i = T_perm[I];
         x[i] /= T.diag[i];
         for (int kk = T.start[i]; kk < T.start[i+1]; kk++){
            x[T.i[kk]] -= T.data[kk] * x[i];
         }
      }
   }
   else {
      for (int I = 0; I < n; I++){
         int i = T_perm[I];
         x[i] = b[i];
         for (int kk = T.start[i]; kk < T.start[i+1]; kk++){
            x[i] -= T.data[kk] * x[T.j[kk]];
         }
         x[i] /= T.diag[i];
      }
   }
}

/****************************
 * Level-scheduled TriSolve 
 ****************************/
void TriSolve_LevelSchedule(TriSolveData *ts,
                            LevelSetData lvl_set, /* level set data */
                            Matrix T, /* triangular matrix */
                            double *x, /* solution (output) */
                            double *b /* right-hand side */
                            )
{
   int lump = 1;

   #pragma omp parallel
   {
      double solve_start;
      int tid = omp_get_thread_num();

      solve_start = omp_get_wtime();

      int num_relax = 0, num_iters = 0;

      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         for (int l = 0; l < lvl_set.num_levels; l++){ /* loop over level sets */
            #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
            for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
               int i = lvl_set.perm[ii];
               x[i] += b[i];
               x[i] /= T.diag[i];
               for (int ii = T.start[i]; ii < T.start[i+1]; ii++){
                  x[T.i[ii]] -= T.data[ii] * x[i];
                  num_relax++;
               }
            }
            num_iters++;
         }
      }
      else {
         for (int l = 0; l < lvl_set.num_levels; l++){ /* loop over level sets */
            #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
            for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
               int i = lvl_set.perm[ii];
               x[i] = b[i];
               for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                  x[i] -= T.data[jj] * x[T.j[jj]];
                  num_relax++;
               }
               x[i] /= T.diag[i];
            }
            num_iters++;
         }
      }
      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;
   }
}

/******************
 * Async TriSolve  
 ******************/
void TriSolve_Async(TriSolveData *ts,
                    Matrix T, /* triangular matrix */
                    double *x, /* solution (output) */
                    double *b /* right-hand side */
                    )
{
   int nnz = T.nnz;
   int n = T.n;
   int lump = 1; 
   int *row_counts;

   int q_size;
   vector<vector<int>> nnz_put_targets(n), col_put_targets(n);
   
   Queue Q;
   /* set up message queue data */
   if (ts->input.MsgQ_flag == 1){
      for (int i = 0; i < n; i++){
         for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
            int j = T.j[jj];
            nnz_put_targets[j].push_back(jj); /* put targets correspond to non-zeros in this row */
            col_put_targets[j].push_back(i);
         }
      }
      q_size = nnz;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
   }
   else {
      /* row counts to track updates to solution vectors */
      row_counts = (int *)calloc(n, sizeof(int));
   }

   #pragma omp parallel
   {
      Matrix T_loc;
      int i_loc, jj_loc, n_loc, nnz_loc, i_prev, i_loc_prev;
      int kk, k;
      double solve_start, setup_start;
      int tid = omp_get_thread_num();
      int *nz_done_flags_loc, *row_done_flags_loc;
      int *row_counts_loc;
      int *my_rows, *my_nzs;
      int num_relax, num_iters;
      int atomic_flag = ts->input.atomic_flag;

      setup_start = omp_get_wtime();

      /* count number of rows and non-zeros assigned to this thread */
      n_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         n_loc++;
      }
      nnz_loc = 0;
      if (ts->input.fine_grained_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int k = 0; k < nnz; k++){
            nnz_loc++;
         }
      }
      else {
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
               nnz_loc++;
            }
         }
      }

      /* for fine-grained solvers, store non-zero indices for this thread*/
      if (ts->input.fine_grained_flag == 1){
         my_nzs = (int *)calloc(nnz_loc+1, sizeof(int));
         kk = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int k = 0; k < nnz; k++){
            my_nzs[kk] = k;
            kk++;
         }
         my_nzs[kk] = nnz;
      }

      /* compute row counts (number of non-zeros per row) and other data.
       * atomics are required for computing row counts in the csc case. */
      i_loc = 0;
      jj_loc = 0;
      if (ts->input.MsgQ_flag == 1){
         row_counts_loc = (int *)calloc(n_loc, sizeof(int));
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            int jj_low = T.start[i];
            int jj_high = T.start[i+1];

            row_counts_loc[i_loc] = jj_high - jj_low;
            i_loc++;
         }
      }
      else {
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
            #pragma omp for schedule(static, lump) nowait
            for (int k = 0; k < nnz; k++){
               #pragma omp atomic
               row_counts[T.i[k]]++;
            }
         }
         else {
            nz_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
            if (ts->input.fine_grained_flag == 0){
               row_done_flags_loc = (int *)calloc(n_loc, sizeof(int));
            }
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               int jj_low = T.start[i];
               int jj_high = T.start[i+1];

               row_counts[i] = jj_high - jj_low;
               if (row_counts[i] == 0 && ts->input.fine_grained_flag == 0){
                  row_done_flags_loc[i_loc] = 1;
               }

               i_loc++;
            }
         }
      }
    
      ts->output.setup_wtime_vec[tid] = omp_get_wtime() - setup_start;
      solve_start = omp_get_wtime();

      /* initialize solution */
      if (ts->input.fine_grained_flag == 1){
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < n; i++){
            x[i] = b[i] / T.diag[i];
         }
      }
      else {
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < n; i++){
            x[i] = b[i];
         }
      }

      num_relax = 0, num_iters = 0;
      jj_loc = 0;
      i_loc = 0;


      /*********************************
       * message queues implementation
       *********************************/
      if (ts->input.MsgQ_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            int jj_start = T.start[i];
            int jj_end = T.start[i+1];
            int jj_diff = jj_end - jj_start;
            while (row_counts_loc[i_loc] > 0){ /* loop until x[i] has been computed */
               int jj_loc_temp = jj_loc;
               for (int jj = jj_start; jj < jj_end; jj++){
                  if (nz_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                     double xj;
                     /* update x[i] if x[j] is available */
                     if (qGet(&Q, jj, &xj)){  
                        x[i] -= T.data[jj] * xj;

                        row_counts_loc[i_loc]--;

                        num_relax++;

                        nz_done_flags_loc[jj_loc_temp] = 1;
                     }
                  }
                  jj_loc_temp++;
               }
               num_iters++;
            }
            x[i] /= T.diag[i];
            /* send x[i] to rows that need it */
            double xi = x[i];
            for (int j = 0; j < nnz_put_targets[i].size(); j++){
               qPut(&Q, nnz_put_targets[i][j], xi);
            }
            jj_loc += jj_diff;
            i_loc++;
         }
      }
      /**************************
       * atomics implementations
       *************************/
      else {
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){ /* compressed sparse column (CSC) version */
            if (ts->input.fine_grained_flag == 1){ /* fine-grained version */
               kk = 0;
               for (int j = 0; j < n; j++){ /* loop over rows */
                  k = my_nzs[kk];
                  int k_end = T.start[j+1];
                  if (k < k_end){ /* does this thread own non-zeros for this row? */
                     int row_counts_j;
                     /* stay idle until element j of x is ready to be used */
                     do {
                        if (atomic_flag == 1){
                           #pragma omp atomic read
                           row_counts_j = row_counts[j];
                        }
                        else {
                           row_counts_j = row_counts[j];
                        }
                     } while (row_counts_j > 0);
                     double xj = x[j];
                     /* for row j, update elements of x and row_counts */
                     while (k < k_end){
                        int i = T.i[k];
                        double Tij = T.data[k];
                        double Tii = T.diag[i];

                        if (atomic_flag == 1){
                           #pragma omp atomic
                           x[i] -= xj * (Tij / Tii);
                        }
                        else {
                           x[i] -= xj * (Tij / Tii);
                        }

                        #pragma omp atomic
                        row_counts[i]--;

                        kk++;
                        k = my_nzs[kk];
                       
                        num_relax++;
                     }
                  }
               }
            }
            else {
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j;
                  /* stay idle until element j of x is ready to be used */
                  do {
                     if (atomic_flag == 1){
                        #pragma omp atomic read
                        row_counts_j = row_counts[j];
                     }
                     else {
                        row_counts_j = row_counts[j];
                     }
                  } while (row_counts_j > 0);
                  /* for row j, update elements of x and row_counts */
                  x[j] /= T.diag[j];
                  double xj = x[j];
                  for (int kk = T.start[j]; kk < T.start[j+1]; kk++){
                     int i = T.i[kk];
                     if (atomic_flag == 1){
                        #pragma omp atomic
                        x[i] -= T.data[kk] * xj;
                     }
                     else {
                        x[i] -= T.data[kk] * xj;
                     }

                     #pragma omp atomic
                     row_counts[i]--;

                     num_relax++;
                  }
               }
            }
         }
         else { /* compressed sparse column (CSR) version */
            if (ts->input.fine_grained_flag == 1){ /* fine-grained version */
               kk = 0;
               for (int i = 0; i < n; i++){ /* loop over rows */
                  k = my_nzs[kk];
                  int kk_start = kk;
                  int done_flag;
                  int k_end = T.start[i+1];
                  if (k < k_end){ /* does this thread own non-zeros for this row? */
                     double Tii = T.diag[i];
                     do{
                        kk = kk_start;
                        k = my_nzs[kk];
                        done_flag = 1;
                        while (k < k_end){
                           if (nz_done_flags_loc[kk] == 0){ /* has this non-zero been used? */
                              int j = T.j[k];
                              int row_counts_j;
                              /* check if x[j] is available (row_counts[j] must be zero) */
                              if (atomic_flag == 1){
                                 #pragma omp atomic read
                                 row_counts_j = row_counts[j];
                              }
                              else {
                                 row_counts_j = row_counts[j];
                              }
                              
                              /* if x[j] is available, update x[i] and row_counts[i] */
                              if (row_counts_j == 0){
                                 double Tij = T.data[k];
                                 double xj = x[j];
                               
                                 if (atomic_flag == 1){
                                    #pragma omp atomic
                                    x[i] -= xj * Tij / Tii;
                                 }
                                 else {
                                    x[i] -= xj * Tij / Tii;
                                 }                           

                                 num_relax++;

                                 #pragma omp atomic
                                 row_counts[i]--;

                                 nz_done_flags_loc[kk] = 1;
                              }
                              else {
                                 done_flag = 0;
                              }
                           }
                           kk++;
                           k = my_nzs[kk];
                        }
                     } while (done_flag == 0); /* loop until x[i] has been completed */
                  }
               }
            } 
            else {
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < n; i++){ /* loop over rows */
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  while (row_done_flags_loc[i_loc] == 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        if (nz_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                           int j = T.j[jj];
                           int row_counts_j;
                           /* check if x[j] is available (row_counts[j] must be zero) */
                           if (atomic_flag == 1){
                              #pragma omp atomic read
                              row_counts_j = row_counts[j];
                           }
                           else {
                              row_counts_j = row_counts[j];
                           }
                           /* if x[j] is available, update x[i] and row_counts[i] */
                           if (row_counts_j == -1){
                              x[i] -= T.data[jj] * x[j];

                              num_relax++;

                              #pragma omp atomic
                              row_counts[i]--;
                              if (row_counts[i] == 0){
                                 row_done_flags_loc[i_loc] = 1;
                              }

                              nz_done_flags_loc[jj_loc_temp] = 1;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  x[i] /= T.diag[i];
                  #pragma omp atomic
                  row_counts[i]--;

                  jj_loc += jj_diff;
                  i_loc++;
               }
            }
         }

         ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

         ts->output.num_relax[tid] = num_relax;
         ts->output.num_iters[tid] = num_iters;

         if (ts->input.MsgQ_flag == 1){

         }
         else {
            if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
               if (ts->input.fine_grained_flag == 1){
                  free(my_nzs);
               }
            }
            else {
               free(nz_done_flags_loc);
               if (ts->input.fine_grained_flag == 1){
                  free(my_nzs);
               }
               else {
                  free(row_done_flags_loc);
               }
            }
         }
      }
   }

   if (ts->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
   else {
      free(row_counts);
   }
}
