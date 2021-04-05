#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/* ******************************************************************************
 * Serial TriSolve (sparse matrix must be in compressed sparse row (CSR) format) 
 * ******************************************************************************/
void TriSolve_CSR(TriSolveData *ts,
                  CSR T, /* triangular matrix */
                  int *T_perm, /* ordering for computing elements of x */
                  double *x, /* solution (output) */
                  double *b /* right-hand side */
                  )
{
   int n = T.n;

   for (int ii = 0; ii < n; ii++){
      int i = T_perm[ii];
      x[i] = b[i] / T.diag[i];
      for (int jj = T.i_ptr[i]; jj < T.i_ptr[i+1]; jj++){
         x[i] -= T.data[jj] * x[T.j[jj]] / T.diag[i];
      }
   }
}

/* ****************************************************************************************
 * Level-scheduled TriSolve (sparse matrix must be in compressed sparse row (CSR) format) 
 * ****************************************************************************************/
void TriSolve_LevelSets_CSR(TriSolveData *ts,
                            LevelSetData lvl_set, /* level set data */
                            CSR T, /* triangular matrix */
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

      for (int l = 0; l < lvl_set.num_levels; l++){ /* loop over level sets */
         #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
         for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
            int i = lvl_set.perm[ii];
            x[i] = b[i] / T.diag[i];
            for (int jj = T.i_ptr[i]; jj < T.i_ptr[i+1]; jj++){
               x[i] -= T.data[jj] * x[T.j[jj]] / T.diag[i];
               num_relax++;
            }
         }
         num_iters++;
      }
      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;
   }
}

/********************************************************************************
 * Async TriSolve (sparse matrix must be in compressed sparse row (CSR) format) 
 * ******************************************************************************/
void TriSolve_Async_CSR(TriSolveData *ts,
                        CSR T, /* triangular matrix */
                        double *x, /* solution (output) */
                        double *b /* right-hand side */
                        )
{
   int nnz = T.nnz;
   int n = T.n;

   int lump = 1;
 
   int *row_counts;

   int q_size;
   vector<vector<int>> put_targets(n);
   Queue Q;
   /* set up message queues */
   if (ts->input.MsgQ_flag == 1){
      for (int i = 0; i < n; i++){
         for (int jj = T.i_ptr[i]; jj < T.i_ptr[i+1]; jj++){
            int ii = T.j[jj];
            put_targets[ii].push_back(jj); /* put targets correspond to non-zeros in this row */
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
      CSR T_loc;
      int i_loc, jj_loc, n_loc, nnz_loc;
      int *T_i_ptr_low, *T_i_ptr_high, *T_i_ptr_diff;
      int accum_count;
      double solve_start, setup_start;
      int tid = omp_get_thread_num();
      int *col_done_flags_loc, *row_done_flags_loc;
      int *row_counts_loc;
      int *my_rows;
      int num_relax, num_iters;

      setup_start = omp_get_wtime();

      /* count number of rows and non-zeros assigned to this thread */
      n_loc = 0;
      nnz_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         n_loc++;
         for (int jj = T.i_ptr[i]; jj < T.i_ptr[i+1]; jj++){
            nnz_loc++;
         }
      }
      /* accum_count is the total number of updates done by this thread.
       * when accum_count is decrimented to zero, convergence is achieved. */
      accum_count = nnz_loc;

      //int block_size = std::max(1, std::min(ts->input.block_size, n_loc));

      col_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));

      //T_i_ptr_low = (int *)calloc(n_loc, sizeof(int));
      //T_i_ptr_high = (int *)calloc(n_loc, sizeof(int));
      //my_rows = (int *)calloc(n_loc, sizeof(int));

      //T_i_ptr_diff = (int *)calloc(n_loc, sizeof(int));

      /* compute row counts (number of non-zeros per row) and other data */
      i_loc = 0;
      jj_loc = 0;
      if (ts->input.MsgQ_flag == 1){
         row_counts_loc = (int *)calloc(n_loc, sizeof(int));

         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            //T_i_ptr_low[i_loc] = T.i_ptr[i];
            //T_i_ptr_high[i_loc] = T.i_ptr[i+1];
            //my_rows[i_loc] = i;
            //T_i_ptr_diff[i_loc] = T_i_ptr_high[i_loc] - T_i_ptr_low[i_loc];

            //T_i_ptr_diff[i_loc] = T.i_ptr[i+1] - T.i_ptr[i];
            //row_counts[i] = T_i_ptr_diff[i_loc];

            int jj_low = T.i_ptr[i];
            int jj_high = T.i_ptr[i+1];

            row_counts_loc[i_loc] = jj_high - jj_low;
            i_loc++;
         }
      }
      else {
         row_done_flags_loc = (int *)calloc(n_loc, sizeof(int));
         T_loc.j = (int *)calloc(nnz_loc, sizeof(int));

         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            //T_i_ptr_low[i_loc] = T.i_ptr[i];
            //T_i_ptr_high[i_loc] = T.i_ptr[i+1];
            //my_rows[i_loc] = i;
            //T_i_ptr_diff[i_loc] = T_i_ptr_high[i_loc] - T_i_ptr_low[i_loc];

            //T_i_ptr_diff[i_loc] = T.i_ptr[i+1] - T.i_ptr[i];
            //row_counts[i] = T_i_ptr_diff[i_loc];

            int jj_low = T.i_ptr[i];
            int jj_high = T.i_ptr[i+1];

            for (int jj = jj_low; jj < jj_high; jj++){
               T_loc.j[jj_loc] = T.j[jj];
               jj_loc++;
            }

            row_counts[i] = jj_high - jj_low;
            if (row_counts[i] == 0){
               row_done_flags_loc[i_loc] = 1;
            }
            i_loc++;
         }
      }
    
      ts->output.setup_wtime_vec[tid] = omp_get_wtime() - setup_start;
      solve_start = omp_get_wtime();

      /* initialize solution */ 
      #pragma omp for schedule(static, lump)
      for (int i = 0; i < n; i++){
         x[i] = b[i] / T.diag[i];
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
            int jj_start = T.i_ptr[i];
            int jj_end = T.i_ptr[i+1];
            int jj_diff = jj_end - jj_start;
            while (row_counts_loc[i_loc] > 0){ /* loop until x[i] has been computed */
               int jj_loc_temp = jj_loc;
               for (int jj = jj_start; jj < jj_end; jj++){
                  if (col_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                     double xj;
                     /* update x[i] if x[j] is available */
                     if (qGet(&Q, jj, &xj)){  
                        x[i] -= T.data[jj] * xj / T.diag[i];

                        row_counts_loc[i_loc]--;

                        accum_count--;

                        num_relax++;

                        col_done_flags_loc[jj_loc_temp] = 1;
                     }
                  }
                  jj_loc_temp++;
               }
            }
            /* send x[i] to rows that need it */
            double xi = x[i];
            for (int j = 0; j < put_targets[i].size(); j++){
               qPut(&Q, put_targets[i][j], xi);
            }
            jj_loc += jj_diff;
            i_loc++;
         }
      }
      /*************************
       * atomics implementation 
       *************************/
      else {
         if (ts->input.atomic_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               int jj_start = T.i_ptr[i];
               int jj_end = T.i_ptr[i+1];
               int jj_diff = jj_end - jj_start;
               while (row_done_flags_loc[i_loc] == 0){ /* loop until x[i] has been completed */
                  int jj_loc_temp = jj_loc;
                  for (int jj = jj_start; jj < jj_end; jj++){
                     if (col_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                        //int j = T.j[jj];
                        int j = T_loc.j[jj_loc_temp];
                        int row_counts_j;
                        /* check if x[j] is available */
                        #pragma omp atomic read
                        row_counts_j = row_counts[j];
                        if (row_counts_j == 0){
                           /* update x and row counts for i */
                           x[i] -= T.data[jj] * x[j] / T.diag[i];

                           #pragma omp atomic
                           row_counts[i]--;
                           if (row_counts[i] == 0){
                              row_done_flags_loc[i_loc] = 1;
                           }

                           accum_count--;

                           num_relax++;

                           col_done_flags_loc[jj_loc_temp] = 1;
                        }
                     }
                     jj_loc_temp++;
                  }
               }
               jj_loc += jj_diff;
               i_loc++;
            }
         }
         else {
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               int jj_start = T.i_ptr[i];
               int jj_end = T.i_ptr[i+1];
               int jj_diff = jj_end - jj_start;
               while (row_done_flags_loc[i_loc] == 0){ /* loop until x[i] has been completed */
                  int jj_loc_temp = jj_loc;
                  for (int jj = jj_start; jj < jj_end; jj++){
                     if (col_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                        //int j = T.j[jj];
                        int j = T_loc.j[jj_loc_temp];
                        int row_counts_j;
                        /* check if x[j] is available */
                        //#pragma omp atomic read
                        row_counts_j = row_counts[j];
                        if (row_counts_j == 0){
                           /* update x and row counts for i */
                           x[i] -= T.data[jj] * x[j] / T.diag[i];

                           #pragma omp atomic
                           row_counts[i]--;
                           if (row_counts[i] == 0){
                              row_done_flags_loc[i_loc] = 1;
                           }

                           accum_count--;

                           num_relax++;

                           col_done_flags_loc[jj_loc_temp] = 1;
                        }
                     }
                     jj_loc_temp++;
                  }
               }
               jj_loc += jj_diff;
               i_loc++;
            }
         }
      }

      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;

      //free(T_i_ptr_low);
      //free(T_i_ptr_high);
      //free(my_rows);
      //free(T_i_ptr_diff);

      if (ts->input.MsgQ_flag == 1){

      }
      else {
         free(T_loc.j);
         free(row_done_flags_loc);
      }

      free(col_done_flags_loc);
   }

   if (ts->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
   else {
      free(row_counts);
   }
}
