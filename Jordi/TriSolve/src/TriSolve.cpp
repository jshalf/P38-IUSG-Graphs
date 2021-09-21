#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

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
   int n = T.n;

   #pragma omp parallel
   {
      double solve_start, solve_stop;
      int tid = omp_get_thread_num();
      int *n_loc, **my_rows, i_loc;

      n_loc = (int *)calloc(lvl_set.num_levels, sizeof(int));
      for (int l = 0; l < lvl_set.num_levels; l++){
         #pragma omp for schedule(static, lump) nowait
         for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
            n_loc[l]++;
         }
      }
      my_rows = (int **)calloc(lvl_set.num_levels, sizeof(int *));
      for (int l = 0; l < lvl_set.num_levels; l++){
         my_rows[l] = (int *)calloc(n_loc[l], sizeof(int *));
         i_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
            int i = lvl_set.perm[ii];
            my_rows[l][i_loc] = i;
            i_loc++;
         }
      }

      int num_relax = 0, num_iters = 0;
      solve_start = omp_get_wtime();

      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < n; i++) x[i] = b[i];
         for (int l = 0; l < lvl_set.num_levels; l++){ /* loop over level sets */
            #pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
            for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
               int i = lvl_set.perm[ii];
               x[i] /= T.diag[i];
               for (int ii = T.start[i]; ii < T.start[i+1]; ii++){
                  #pragma omp atomic
                  x[T.i[ii]] -= T.data[ii] * x[i];
                  num_relax++;
               }
            }
            num_iters++;
         }
      }
      else {
         i_loc = 0;
         for (int l = 0; l < lvl_set.num_levels; l++){ /* loop over level sets */
            //#pragma omp for schedule(static, lump) /* parallel loop over elements within level set */
            //for (int ii = lvl_set.level_start[l]; ii < lvl_set.level_start[l+1]; ii++){
            for (int i_loc = 0; i_loc < n_loc[l]; i_loc++){
               //int i = lvl_set.perm[ii];
               int i = my_rows[l][i_loc];
               x[i] = b[i];
               for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                  x[i] -= T.data[jj] * x[T.j[jj]];
                  num_relax++;
               }
               x[i] /= T.diag[i];
            }
            num_iters++;
            #pragma omp barrier
            num_iters++;
         }
      }
      solve_stop = omp_get_wtime();
      ts->output.solve_wtime_vec[tid] = solve_stop - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;

      for (int l = 0; l < lvl_set.num_levels; l++){
         free(my_rows[l]);
      }
      free(my_rows);
      free(n_loc);
   }
}

void TriSolve_AtomicCounter(TriSolveData *ts, Matrix T, /* triangular matrix */
                            double *x,                  /* solution (output) */
                            double *b                   /* right-hand side */
)
{
   int n = T.n;
   cache *counters = (cache *)calloc(n, sizeof(cache));

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      double solve_start = omp_get_wtime();
      #pragma omp for schedule(static, 1)
      for (int i = 0; i < n; i++) {
         double z = b[i];
         for (int kk = T.start[i]; kk < T.start[i + 1]; kk++) {
            int idx = T.j[kk];
            int rdy;
            do {
               #pragma omp atomic read
               rdy = counters[idx].c;
            } while (!rdy);
            z += -1 * (T.data[kk] * x[idx]);
         }
         x[i] = z / T.diag[i];
         #pragma omp atomic write
         counters[i].c = 1;
      }
      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;
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
   int *row_counts, *row_done_flags, **row_counts_thread;
   LevelSetData lvl_set = ts->L_lvl_set;
   //cache *row_done_flags;
   
   Queue Q, Q_comm;
   vector<int> row_to_thread;
   int *t_sum, **t_sum_thread;
   if (ts->input.MsgQ_flag == 1){ /* set up message queue data */
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         row_counts = (int *)calloc(n, sizeof(int));
         row_counts_thread = (int **)calloc(ts->input.num_threads, sizeof(int *));
         for (int t = 0; t < ts->input.num_threads; t++){
            row_counts_thread[t] = (int *)calloc(n, sizeof(int));
         }

         qAlloc(&Q, n);
      }
      else {
         row_to_thread.resize(n);
         t_sum = (int *)calloc(ts->input.num_threads, sizeof(int));
         t_sum_thread = (int **)calloc(ts->input.num_threads, sizeof(int *));
         for (int t = 0; t < ts->input.num_threads; t++){
            t_sum_thread[t] = (int *)calloc(ts->input.num_threads, sizeof(int *));
         }

         qAlloc(&Q_comm, n);
         qInitLock(&Q_comm);

         qAlloc(&Q, nnz);
      }

      qInitLock(&Q);
   }
   else {
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         row_counts = (int *)calloc(n, sizeof(int));
      }
      else {
         if (ts->input.fine_grained_flag == 1){
            row_counts = (int *)calloc(n, sizeof(int));
         }
         else {
            //row_done_flags = (cache *)calloc(n, sizeof(cache));
            row_done_flags = (int *)calloc(n, sizeof(int));
         }
      }
   }

   int *row_perm, *nz_perm;
   if (ts->input.fine_grained_flag == 1){
      row_perm = (int *)calloc(n, sizeof(int));
      nz_perm = (int *)calloc(nnz, sizeof(int));
      if (ts->input.solver_type == TRISOLVE_ASYNC_LEVEL_SCHEDULED){
         int k = 0; 
         for (int ii = 0; ii < n; ii++){
            int i = lvl_set.perm[ii];
            for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
               nz_perm[k] = jj;
               k++;
            }
         }
      }
   }

   #pragma omp parallel
   {
      Matrix T_loc;
      int i_loc, j_loc, jj_loc, n_loc, nnz_loc, i_prev, i_loc_prev;
      int kk, k;
      double solve_wtime, setup_wtime, wtime_start, wtime_stop;
      double comp_wtime_start, comp_wtime_stop, comp_wtime = 0.0;
      double MsgQ_wtime_start, MsgQ_wtime_stop, MsgQ_wtime = 0.0;
      double MsgQ_put_wtime_start, MsgQ_put_wtime_stop, MsgQ_put_wtime = 0.0;
      double MsgQ_get_wtime_start, MsgQ_get_wtime_stop, MsgQ_get_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles_stop, MsgQ_cycles = 0;
      uint64_t MsgQ_put_cycles_start, MsgQ_put_cycles_stop, MsgQ_put_cycles = 0;
      uint64_t MsgQ_get_cycles_start, MsgQ_get_cycles_stop, MsgQ_get_cycles = 0;
      int num_qGets = 0, num_qPuts = 0;
      int dummy_num_qGets = 0, dummy_num_qPuts = 0;
      int dummy_num_spins = 0;
      int temp_num_qPuts = 0, temp_num_qGets = 0;
      int tid = omp_get_thread_num();
      int *nz_done_flags_loc, *row_done_flags_loc;
      int *row_counts_loc;
      int *my_rows, *my_nzs;
      int num_relax, num_iters;
      int atomic_flag = ts->input.atomic_flag;
      double *x_loc;
      double dummy;
      int num_gets;
      vector<vector<int>> put_targets;

      wtime_start = omp_get_wtime();

      n_loc = 0;
      nnz_loc = 0;
      if (ts->input.fine_grained_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int k = 0; k < nnz; k++){
            nnz_loc++;
         }
         my_nzs = (int *)calloc(nnz_loc+1, sizeof(int));
         if (ts->input.solver_type == TRISOLVE_ASYNC_LEVEL_SCHEDULED){
            kk = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int k = 0; k < nnz; k++){
               my_nzs[kk] = nz_perm[k];
               kk++;
            }
            my_nzs[kk] = nz_perm[nnz-1];
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               row_perm[i] = lvl_set.perm[i];
            }
         }
         else {
            kk = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int k = 0; k < nnz; k++){
               my_nzs[kk] = k;
               kk++;
            }
            my_nzs[kk] = nnz;
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               row_perm[i] = i;
            }
         }
      }
      else {
         if (ts->input.solver_type == TRISOLVE_ASYNC_LEVEL_SCHEDULED){ 
            #pragma omp for schedule(static, lump) nowait
            for (int ii = 0; ii < n; ii++){
               n_loc++;
               int i = lvl_set.perm[ii];
               for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                  nnz_loc++;
               }
            }
            my_rows = (int *)calloc(n_loc, sizeof(int));
            i_loc = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int ii = 0; ii < n; ii++){
               my_rows[i_loc] = lvl_set.perm[ii];
               i_loc++;
            }
         }
         else {
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               n_loc++;
               for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                  nnz_loc++;
               }
            }
            my_rows = (int *)calloc(n_loc, sizeof(int));
            i_loc = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               my_rows[i_loc] = i;
               i_loc++;
            }
         }

         x_loc = (double *)calloc(n_loc, sizeof(double));
         T_loc.diag = (double *)calloc(n_loc, sizeof(double));
         T_loc.start = (int *)calloc(n_loc+1, sizeof(int));
         T_loc.data = (double *)calloc(nnz_loc, sizeof(double));
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
            T_loc.i = (int *)calloc(nnz_loc, sizeof(int));
         }
         else {
            T_loc.j = (int *)calloc(nnz_loc, sizeof(int));
         }
         i_loc = 0, jj_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            T_loc.diag[i_loc] = T.diag[i];
            int i_nnz = 0;
            for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
               if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
                  T_loc.i[jj_loc] = T.i[jj];
               }
               else {
                  T_loc.j[jj_loc] = T.j[jj];
               }
               T_loc.data[jj_loc] = T.data[jj];
               jj_loc++;
               i_nnz++;
            }
            T_loc.start[i_loc+1] = T_loc.start[i_loc] + i_nnz;
            i_loc++;
         }
      }

      /* compute row counts (number of non-zeros per row) and other data.
       * atomics are required for computing row counts in the csc case. */
      i_loc = 0;
      jj_loc = 0;
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         if (ts->input.MsgQ_flag == 1){
            #pragma omp for schedule(static, lump)
            for (int k = 0; k < nnz; k++){
               row_counts_thread[tid][T.i[k]]++;
            }
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               for (int t = 0; t < ts->input.num_threads; t++){
                  row_counts[i] += row_counts_thread[t][i];
               }
            }
         }
         else {
            #pragma omp for schedule(static, lump) nowait
            for (int k = 0; k < nnz; k++){
               #pragma omp atomic
               row_counts[T.i[k]]++;
            }
         }
      }
      else {
         if (ts->input.fine_grained_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               int jj_low = T.start[i];
               int jj_high = T.start[i+1];

               row_counts[i] = jj_high - jj_low;
            }
         }
         else {

            if (ts->input.MsgQ_flag == 1){
               put_targets.resize(n_loc);

               for (i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  row_to_thread[i] = tid;
               }
               #pragma omp barrier
               for (i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                     int ii = T.j[jj];
                     t_sum_thread[tid][row_to_thread[ii]]++;
                  }
               }
               #pragma omp barrier
               #pragma omp for schedule(static, lump)
               for (int t = 0; t < ts->input.num_threads; t++){
                  for (int tt = 0; tt < ts->input.num_threads; tt++){
                     t_sum[t] += t_sum_thread[tt][t];
                  }
               }
               int num_get = t_sum[tid];
               int temp_count = 0;
               for (i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  for (int jj = T.start[i]; jj < T.start[i+1]; jj++){
                     int ii = T.j[jj];
                     qPut(&Q_comm, ii, (double)jj);
                     temp_count++;
                  }
               }
               temp_count = 0;
               while (temp_count < num_get) {
                  for (i_loc = 0; i_loc < n_loc; i_loc++){
                     int i = my_rows[i_loc];
                     double target;
                     if (qGet(&Q_comm, i, &target)){
                        put_targets[i_loc].push_back((int)target);
                        temp_count++;
                     }
                  }
               }
               num_gets = 0;
               for (i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  num_gets += jj_end - jj_start;
               }
            }
            else {
            }
         }
         nz_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
      }

      /* initialize solution */
      if (ts->input.MsgQ_flag == 1){
         if (ts->input.fine_grained_flag == 1){
            i_loc = 0;
            #pragma omp for schedule(static, lump)
            for (int i = 0; i < n; i++){
               x_loc[i_loc] = b[i] / T.diag[i];
               i_loc++;
            }
         }
         else {
            for (i_loc = 0; i_loc < n_loc; i_loc++){
               int i = my_rows[i_loc];
               x_loc[i_loc] = b[i];
            }
            #pragma omp barrier
         }
      }
      else {
         if (ts->input.fine_grained_flag == 1){
            #pragma omp for schedule(static, lump)
            for (int i = 0; i < n; i++){
               x[i] = b[i] / T.diag[i];
            }
         }
         else {
            for (i_loc = 0; i_loc < n_loc; i_loc++){
               int i = my_rows[i_loc];
               x[i] = b[i];
            }
            #pragma omp barrier
         }
      }

      num_relax = 0, num_iters = 0;
      jj_loc = 0;
      i_loc = 0;

      wtime_stop = omp_get_wtime();
      ts->output.setup_wtime_vec[tid] = wtime_stop - wtime_start;








      wtime_start = omp_get_wtime();
/********************************************************************
 *            
 *                       MESSAGE QUEUES
 *
 ********************************************************************/
      if (ts->input.MsgQ_flag == 1){

         /**********************
          *
          *      MsgQ CSC
          *
          **********************/         
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
            #pragma omp for schedule(static, lump) nowait
            for (int j = 0; j < n; j++){ /* loop over rows */
               int row_counts_j = row_counts[j];
               double z;
               /* stay idle until element j of x is ready to be used */
               while (row_counts_j > 0){
                  int get_flag = qGet(&Q, j, &z);
                  if (get_flag == 1){
                     x_loc[j_loc] -= z;
                     row_counts_j--;
                  }
               }
               /* for row j, update elements of x and row_counts */
               x_loc[j_loc] /= T.diag[j];
               double xj = x_loc[j_loc];
               for (int kk = T.start[j]; kk < T.start[j+1]; kk++){
                  int i = T.i[kk];
                  z = T.data[kk] * xj;
                  qPut(&Q, i, z);
                  num_relax++;
               }
               j_loc++;
            }
         }










         /**********************
          *
          *      MsgQ CSR
          *
          **********************/
         else {
            //jj_loc = 0;
            //for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
            //   int i = my_rows[i_loc];
            //   int jj_start = T.start[i];
            //   int jj_end = T.start[i+1];
            //   int jj_diff = jj_end - jj_start;
            //   int row_count_i = jj_diff;
            //   while (row_count_i > 0){ /* loop until x[i] has been completed */
            //      int jj_loc_temp = jj_loc;
            //      for (int jj = jj_start; jj < jj_end; jj++){
            //         int this_nz_done = nz_done_flags_loc[jj_loc_temp];
            //         if (this_nz_done == 0){ /* has this non-zero been used? */
            //            double xj;
            //            int get_flag = qGet(&Q, jj, &xj);
            //            /* if x[j] is available, update x[i] and other data */
            //            if (get_flag == 1){
            //               x_loc[i_loc] -= T_loc.data[jj_loc_temp] * xj;
            //               row_count_i--;
            //               nz_done_flags_loc[jj_loc_temp] = 1;
            //               num_qGets++;
            //            }
            //         }
            //         jj_loc_temp++;
            //      }
            //   }
            //   x_loc[i_loc] /= T_loc.diag[i_loc];
            //   double xi = x_loc[i_loc];
            //   int num_puts_i = put_targets[i_loc].size();
            //   for (int j = 0; j < num_puts_i; j++){
            //      int put_target = put_targets[i_loc][j];
            //      qPut(&Q, put_target, xi);
            //      num_qPuts++;
            //   }

            //   jj_loc += jj_diff;
            //}
            //for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
            //   nz_done_flags_loc[jj_loc] = 0;
            //}
            //#pragma omp barrier

            /*****************
             *   MsgQ wtime
             *****************/
            if (ts->input.MsgQ_wtime_flag == 1){
               jj_loc = 0;
               MsgQ_wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj;
                           int get_flag = qGet(&Q, jj, &xj);
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += xj;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  double xi =1.0;
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     qPut(&Q, put_target, xi);
                     num_qPuts++;
                  }

                  jj_loc += jj_diff;
               }
               MsgQ_wtime_stop = omp_get_wtime();
               MsgQ_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;

               for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
                  nz_done_flags_loc[jj_loc] = 0;
               }
               #pragma omp barrier

               jj_loc = 0;
               MsgQ_wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj = 1.0;
                           int get_flag = 1;
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += (double)get_flag + (double)jj + xj;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  double xi =1.0;
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     num_qPuts++;
                     dummy += (double)put_target;
                  }

                  jj_loc += jj_diff;
               }
               MsgQ_wtime_stop = omp_get_wtime();
               MsgQ_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;

               MsgQ_put_wtime = MsgQ_wtime;
            }
            /*****************
             *   MsgQ cycles
             *****************/
            else if (ts->input.MsgQ_cycles_flag == 1){
            }
           /*****************
            *   comp wtime
            *****************/
            else if (ts->input.comp_wtime_flag == 1){
               jj_loc = 0;
               comp_wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj = 1.0;
                           int get_flag = 1;
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              x_loc[i_loc] -= T_loc.data[jj_loc_temp] * xj;
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += (double)get_flag + (double)jj + xj;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  x_loc[i_loc] /= T_loc.diag[i_loc];
                  double xi = x_loc[i_loc];
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     num_qPuts++;
                     dummy += (double)put_target + xi;
                  }

                  jj_loc += jj_diff;
               }
               comp_wtime_stop = omp_get_wtime();
               comp_wtime += comp_wtime_stop - comp_wtime_start;
            }
           /******************************
            *   MsgQ no-op and comp no-op
            ******************************/
            else if ((ts->input.MsgQ_noop_flag == 1) && (ts->input.comp_noop_flag == 1)){
               jj_loc = 0;
               wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj = 1.0;
                           int get_flag = 1;
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += (double)(xj + get_flag);
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  double xi = 1.0;
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     num_qPuts++;
                     dummy += (double)put_target;
                  }

                  jj_loc += jj_diff;
               }
            }
           /******************
            *   MsgQ no-op
            ******************/
            else if (ts->input.MsgQ_noop_flag == 1){
               jj_loc = 0;
               wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj;
                           int get_flag = 1;
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              x_loc[i_loc] -= T_loc.data[jj_loc_temp] * xj;
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += (double)get_flag;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  x_loc[i_loc] /= T_loc.diag[i_loc];
                  double xi = x_loc[i_loc];
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     num_qPuts++;
                     dummy += (double)put_target + xi;
                  }

                  jj_loc += jj_diff;
               }
            }
           /*****************
            *   comp no-op
            *****************/
            else if (ts->input.comp_noop_flag == 1){
               jj_loc = 0;
               wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj;
                           int get_flag = qGet(&Q, jj, &xj);
                           if (get_flag == 1){
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                              dummy += (double)get_flag + xj;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  double xi = 1.0;
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     qPut(&Q, put_target, xi);
                     num_qPuts++;
                  }

                  jj_loc += jj_diff;
               }
            }
            /**********************************************
             * standard scheme (no timers, no-ops, etc...)
             **********************************************/
            else {
               jj_loc = 0;
               wtime_start = omp_get_wtime();
               for (i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                  int i = my_rows[i_loc];
                  int jj_start = T.start[i];
                  int jj_end = T.start[i+1];
                  int jj_diff = jj_end - jj_start;
                  int row_count_i = jj_diff;
                  while (row_count_i > 0){ /* loop until x[i] has been completed */
                     int jj_loc_temp = jj_loc;
                     for (int jj = jj_start; jj < jj_end; jj++){
                        int this_nz_done = nz_done_flags_loc[jj_loc_temp];
                        if (this_nz_done == 0){ /* has this non-zero been used? */
                           double xj;
                           int get_flag = qGet(&Q, jj, &xj);
                           /* if x[j] is available, update x[i] and other data */
                           if (get_flag == 1){
                              x_loc[i_loc] -= T_loc.data[jj_loc_temp] * xj;
                              row_count_i--;
                              nz_done_flags_loc[jj_loc_temp] = 1;
                              num_qGets++;
                           }
                        }
                        jj_loc_temp++;
                     }
                  }
                  x_loc[i_loc] /= T_loc.diag[i_loc];
                  double xi = x_loc[i_loc];
                  int num_puts_i = put_targets[i_loc].size();
                  for (int j = 0; j < num_puts_i; j++){
                     int put_target = put_targets[i_loc][j];
                     qPut(&Q, put_target, xi);
                     num_qPuts++;
                  }

                  jj_loc += jj_diff;
               }
            }
         }
      }




















/*****************************************************************
 * 
 *                         ATOMIC
 *
 *****************************************************************/
      else {
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
            /**************************
             * 
             * Atomic CSC fine-grained
             *
             **************************/ 
            if (ts->input.fine_grained_flag == 1){
               kk = 0;
               /************
                * atomic
                ************/
               if (atomic_flag == 1){
                  for (int J = 0; J < n; J++){ /* loop over rows */
                     int j = row_perm[J];
                     k = my_nzs[kk];
                     int k_start = T.start[j]; 
                     int k_end = T.start[j+1];
                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
                        int row_counts_j;
                        /* stay idle until element j of x is ready to be used */
                        do {
                           #pragma omp atomic read
                           row_counts_j = row_counts[j];
                        } while (row_counts_j > 0);
                        double xj = x[j];
                        /* for row j, update elements of x and row_counts */
                        while (k < k_end && k >= k_start){
                           int i = T.i[k];
                           double Tij = T.data[k];
                           double Tii = T.diag[i];

                           #pragma omp atomic
                           x[i] -= xj * (Tij / Tii);

                           #pragma omp atomic
                           row_counts[i]--;

                           kk++;
                           k = my_nzs[kk];
                          
                           num_relax++;
                        }
                     }
                  }
               }
               /************
                * no atomic
                ************/
               else {
                  for (int J = 0; J < n; J++){ /* loop over rows */
                     int j = row_perm[J];
                     k = my_nzs[kk];
                     int k_start = T.start[j];
                     int k_end = T.start[j+1];
                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
                        int row_counts_j;
                        /* stay idle until element j of x is ready to be used */
                        do {
                           //#pragma omp atomic read
                           row_counts_j = row_counts[j];
                        } while (row_counts_j > 0);
                        double xj = x[j];
                        /* for row j, update elements of x and row_counts */
                        while (k < k_end && k >= k_start){
                           int i = T.i[k];
                           double Tij = T.data[k];
                           double Tii = T.diag[i];

                           x[i] -= xj * (Tij / Tii);

                           #pragma omp atomic
                           row_counts[i]--;

                           kk++;
                           k = my_nzs[kk];

                           num_relax++;
                        }
                     }
                  }
               }
            }










            /**********************
             * 
             *     Atomic CSC
             *
             **********************/
            else {
               /************
                * atomic
                ************/
               if (atomic_flag == 1){
                  #pragma omp for schedule(static, lump) nowait
                  for (int j = 0; j < n; j++){ /* loop over rows */
                     int row_counts_j;
                     /* stay idle until element j of x is ready to be used */
                     do {
                        #pragma omp atomic read
                        row_counts_j = row_counts[j];
                     } while (row_counts_j > 0);
                     /* for row j, update elements of x and row_counts */
                     x[j] /= T.diag[j];
                     double xj = x[j];
                     for (int kk = T.start[j]; kk < T.start[j+1]; kk++){
                        int i = T.i[kk];

                        #pragma omp atomic
                        x[i] -= T.data[kk] * xj;

                        #pragma omp atomic
                        row_counts[i]--;

                        num_relax++;
                     }
                  }
               }
               /************ 
                * no atomic 
                ************/
               else {
                  #pragma omp for schedule(static, lump) nowait
                  for (int j = 0; j < n; j++){ /* loop over rows */
                     int row_counts_j;
                     /* stay idle until element j of x is ready to be used */
                     do {
                        //#pragma omp atomic read
                        row_counts_j = row_counts[j];
                     } while (row_counts_j > 0);
                     /* for row j, update elements of x and row_counts */
                     x[j] /= T.diag[j];
                     double xj = x[j];
                     for (int kk = T.start[j]; kk < T.start[j+1]; kk++){
                        int i = T.i[kk];

                        x[i] -= T.data[kk] * xj;

                        #pragma omp atomic
                        row_counts[i]--;

                        num_relax++;
                     }
                  }
               }
            }


         }
         else {










            /**************************
             * 
             * Atomic CSR fine-grained
             *
             **************************/
            if (ts->input.fine_grained_flag == 1){
               kk = 0;
               /************
                * atomic
                ************/
               if (atomic_flag == 1){
                  for (int I = 0; I < n; I++){ /* loop over rows */
                     int i = row_perm[I];
                     k = my_nzs[kk];
                     int kk_start = kk;
                     int done_flag;
                     int k_start = T.start[i];
                     int k_end = T.start[i+1];
                     if (k < k_end && k >= k_start){ /* does this thread own non-zeros for this row? */
                        double Tii = T.diag[i];
                        do{
                           kk = kk_start;
                           k = my_nzs[kk];
                           done_flag = 1;
                           while (k < k_end && k >= k_start){
                              if (nz_done_flags_loc[kk] == 0){ /* has this non-zero been used? */
                                 int j = T.j[k];
                                 int row_counts_j;
                                 /* check if x[j] is available (row_counts[j] must be zero) */
                                 #pragma omp atomic read
                                 row_counts_j = row_counts[j];
                                 
                                 /* if x[j] is available, update x[i] and row_counts[i] */
                                 if (row_counts_j == 0){
                                    double Tij = T.data[k];
                                    double xj = x[j];
                                  
                                    #pragma omp atomic
                                    x[i] -= xj * Tij / Tii;

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
               /************
                * no atomic
                ************/
               else {
                  for (int I = 0; I < n; I++){ /* loop over rows */
                     int i = row_perm[I];
                     k = my_nzs[kk];
                     int kk_start = kk;
                     int done_flag;
                     int k_start = T.start[i];
                     int k_end = T.start[i+1];
                     if (k < k_end){ /* does this thread own non-zeros for this row? */
                        double Tii = T.diag[i];
                        do{
                           kk = kk_start;
                           k = my_nzs[kk];
                           done_flag = 1;
                           while (k < k_end && k >= k_start){
                              if (nz_done_flags_loc[kk] == 0){ /* has this non-zero been used? */
                                 int j = T.j[k];
                                 int row_counts_j;
                                 /* check if x[j] is available (row_counts[j] must be zero) */
                                 row_counts_j = row_counts[j];

                                 /* if x[j] is available, update x[i] and row_counts[i] */
                                 if (row_counts_j == 0){
                                    double Tij = T.data[k];
                                    double xj = x[j];

                                    x[i] -= xj * Tij / Tii;

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
            } 










            /**********************
             *        
             *    Atomic CSR
             *
             **********************/
            else {
               /************
                * atomic
                ************/
               if (atomic_flag == 1){
                  for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                     int i = my_rows[i_loc];
                     int jj_start = T.start[i];
                     int jj_end = T.start[i+1];
                     int jj_diff = jj_end - jj_start;
                     int row_count_i = jj_diff;
                     while (row_count_i > 0){ /* loop until x[i] has been completed */
                        int jj_loc_temp = jj_loc;
                        for (int jj = jj_start; jj < jj_end; jj++){
                           if (nz_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                              int j = T.j[jj];
                              int row_done_flag_j;
                              /* check if x[j] is available (row_counts[j] must be zero) */
                              #pragma omp atomic read
                              row_done_flag_j = row_done_flags[j];//row_done_flag_j = row_done_flags[j].c;

                              /* if x[j] is available, update x[i] and row_counts[i] */
                              if (row_done_flag_j == 1){
                                 x[i] -= T.data[jj] * x[j];

                                 num_relax++;

                                 row_count_i--;

                                 nz_done_flags_loc[jj_loc_temp] = 1;
                              }
                           }
                           jj_loc_temp++;
                        }
                     }
                     x[i] /= T.diag[i];
                     #pragma omp atomic write
                     row_done_flags[i] = 1;//row_done_flags[i].c = 1;

                     jj_loc += jj_diff;
                  }
               }
               /************
                * no atomic
                ************/
               else {
                  for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
                     int i = my_rows[i_loc];
                     int jj_start = T.start[i];
                     int jj_end = T.start[i+1];
                     int jj_diff = jj_end - jj_start;
                     int row_count_i = jj_diff;
                     while (row_count_i > 0){ /* loop until x[i] has been completed */
                        int jj_loc_temp = jj_loc;
                        for (int jj = jj_start; jj < jj_end; jj++){
                           if (nz_done_flags_loc[jj_loc_temp] == 0){ /* has this non-zero been used? */
                              int j = T.j[jj];
                              int row_done_flag_j;
                              /* check if x[j] is available (row_counts[j] must be zero) */
                              row_done_flag_j = row_done_flags[j];//row_done_flag_j = row_done_flags[j].c;

                              /* if x[j] is available, update x[i] and row_counts[i] */
                              if (row_done_flag_j == 1){
                                 x[i] -= T.data[jj] * x[j];

                                 num_relax++;

                                 row_count_i--;

                                 nz_done_flags_loc[jj_loc_temp] = 1;
                              }
                           }
                           jj_loc_temp++;
                        }
                     }
                     x[i] /= T.diag[i];
                     row_done_flags[i] = 1;//row_done_flags[i].c = 1;

                     jj_loc += jj_diff;
                  }
               }
            }
         }
      }

      wtime_stop = omp_get_wtime();
      ts->output.solve_wtime_vec[tid] = wtime_stop - wtime_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;

      if (ts->input.MsgQ_flag == 1){
         for (i_loc = 0; i_loc < n_loc; i_loc++){
            int i = my_rows[i_loc];
            x[i] = x_loc[i_loc];
         }

         if (ts->input.MsgQ_wtime_flag == 1){
            ts->output.MsgQ_put_wtime_vec[tid] = MsgQ_put_wtime;
            ts->output.MsgQ_get_wtime_vec[tid] = MsgQ_get_wtime;
         }
         else if (ts->input.MsgQ_cycles_flag == 1){
            ts->output.MsgQ_put_cycles_vec[tid] = MsgQ_put_cycles;
            ts->output.MsgQ_get_cycles_vec[tid] = MsgQ_get_cycles;
         }
         else if (ts->input.comp_wtime_flag == 1){
            ts->output.comp_wtime_vec[tid] = comp_wtime;
         }

         PrintDummy(dummy);
      }

      if (ts->input.fine_grained_flag == 1){
         free(my_nzs);
      }
      else {
         free(my_rows);
         free(x_loc);
         free(T_loc.diag);
         free(T_loc.data);
         free(T_loc.start);
         if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
            free(T_loc.i);
         }
         else {
            free(T_loc.j);
            free(nz_done_flags_loc);
         }
      }
   }

   if (ts->input.MsgQ_flag == 1){
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         free(row_counts);
         for (int t = 0; t < ts->input.num_threads; t++){
            free(row_counts_thread[t]);
         }
         free(row_counts_thread);
      }
      else {
         free(t_sum);
         for (int t = 0; t < ts->input.num_threads; t++){
            free(t_sum_thread[t]);
         }
         free(t_sum_thread);
         qDestroyLock(&Q_comm);
         qFree(&Q_comm);
      }
      qDestroyLock(&Q);
      qFree(&Q);
   }
   else {
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         free(row_counts);
      }
      else {
         if (ts->input.fine_grained_flag == 1){
            free(row_counts);
         }
         else {
            free(row_done_flags);
         }
      }
   }

   if (ts->input.fine_grained_flag == 1){
      free(nz_perm);
      free(row_perm);
   }

   //double min_elem = *min_element(ts->output.solve_wtime_vec, ts->output.solve_wtime_vec+ts->input.num_threads);
   //double max_elem = *max_element(ts->output.solve_wtime_vec, ts->output.solve_wtime_vec+ts->input.num_threads);
   //printf("%e %e\n", min_elem, max_elem);
}
