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
   int *row_counts, *row_done_flags;
   //cache *row_done_flags;

   int q_size;
   
   Queue Q;
   /* set up message queue data */
   if (ts->input.MsgQ_flag == 1){
      q_size = nnz;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
   }

   /* row counts to track updates to solution vectors */
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

   #pragma omp parallel
   {
      Matrix T_loc;
      int i_loc, j_loc, jj_loc, n_loc, nnz_loc, i_prev, i_loc_prev;
      int kk, k;
      double solve_start, setup_start;
      double comp_wtime_start, comp_wtime = 0.0, MsgQ_wtime_start, MsgQ_wtime = 0.0; 
      uint64_t MsgQ_cycles_start, MsgQ_cycles = 0, comp_cycles_start, comp_cycles = 0;
      int tid = omp_get_thread_num();
      int *nz_done_flags_loc, *row_done_flags_loc;
      int *row_counts_loc;
      int *my_rows, *my_nzs;
      int num_relax, num_iters;
      int atomic_flag = ts->input.atomic_flag;
      double *x_loc;

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
      if (ts->input.mat_storage_type == MATRIX_STORAGE_CSC){
         #pragma omp for schedule(static, lump) nowait
         for (int k = 0; k < nnz; k++){
            #pragma omp atomic
            row_counts[T.i[k]]++;
         }
      }
      else {
         nz_done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
         if (ts->input.fine_grained_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < n; i++){
               int jj_low = T.start[i];
               int jj_high = T.start[i+1];

               row_counts[i] = jj_high - jj_low;
            }
         }
         else {
         }
      }
    
      ts->output.setup_wtime_vec[tid] = omp_get_wtime() - setup_start;
      solve_start = omp_get_wtime();

      /* initialize solution */
      if (ts->input.MsgQ_flag == 1){
         x_loc = (double *)calloc(n_loc, sizeof(double));
         i_loc = 0;
         if (ts->input.fine_grained_flag == 1){
            #pragma omp for schedule(static, lump)
            for (int i = 0; i < n; i++){
               x_loc[i_loc] = b[i] / T.diag[i];
               i_loc++;
            }
         }
         else {
            #pragma omp for schedule(static, lump)
            for (int i = 0; i < n; i++){
               x_loc[i_loc] = b[i];
               i_loc++;
            }
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
            #pragma omp for schedule(static, lump)
            for (int i = 0; i < n; i++){
               x[i] = b[i];
            }
         }
      }

      num_relax = 0, num_iters = 0;
      jj_loc = 0;
      i_loc = 0;




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
            j_loc = 0;
           /*****************
            *   MsgQ wtime
            *****************/
            if (ts->input.MsgQ_wtime_flag == 1){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z;
                  /* stay idle until element j of x is ready to be used */
                  while (row_counts_j > 0){
                     MsgQ_wtime_start = omp_get_wtime();
                     int get_flag = qGet(&Q, j, &z);
                     MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                     if (get_flag == 1){
                        x_loc[j_loc] -= z;
                        row_counts_j--;
                     }
                  }
                  /* for row j, update elements of x and row_counts */
                  x_loc[j_loc] /= T.diag[j];
                  double xj = x_loc[j_loc];
                  int low = T.start[j];
                  int high = T.start[j+1];
                  for (int kk = low; kk < high; kk++){
                     int i = T.i[kk];
                     z = T.data[kk] * xj;
                     MsgQ_wtime_start = omp_get_wtime();
                     qPut(&Q, i, z);
                     MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                  }
                  j_loc++;
               }
            }
            /*****************
            *   MsgQ cycles
            *****************/
            else if (ts->input.MsgQ_cycles_flag == 1){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z;
                  /* stay idle until element j of x is ready to be used */
                  while (row_counts_j > 0){
                     MsgQ_cycles_start = rdtsc();
                     int get_flag = qGet(&Q, j, &z);
                     MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                     if (get_flag == 1){
                        x_loc[j_loc] -= z;
                        row_counts_j--;
                     }
                  }
                  /* for row j, update elements of x and row_counts */
                  x_loc[j_loc] /= T.diag[j];
                  double xj = x_loc[j_loc];
                  int low = T.start[j];
                  int high = T.start[j+1];
                  for (int kk = low; kk < high; kk++){
                     int i = T.i[kk];
                     z = T.data[kk] * xj;
                     MsgQ_cycles_start = rdtsc();
                     qPut(&Q, i, z);
                     MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                  }
                  j_loc++;
               }
            }
           /*****************
            *   comp wtime
            *****************/
            else if (ts->input.comp_wtime_flag == 1){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z;
                  /* stay idle until element j of x is ready to be used */
                  while (row_counts_j > 0){
                     int get_flag = qGet(&Q, j, &z);
                     if (get_flag == 1){
                        comp_wtime_start = omp_get_wtime();
                        x_loc[j_loc] -= z;
                        comp_wtime += omp_get_wtime() - comp_wtime_start;
                        row_counts_j--;
                     }
                  }
                  comp_wtime_start = omp_get_wtime();
                  /* for row j, update elements of x and row_counts */
                  x_loc[j_loc] /= T.diag[j];
                  double xj = x_loc[j_loc];
                  int low = T.start[j];
                  int high = T.start[j+1];
                  comp_wtime += omp_get_wtime() - comp_wtime_start;
                  for (int kk = low; kk < high; kk++){
                     comp_wtime_start = omp_get_wtime();
                     int i = T.i[kk];
                     z = T.data[kk] * xj;
                     comp_wtime += omp_get_wtime() - comp_wtime_start;
                     qPut(&Q, i, z);
                  }
                  j_loc++;
               }
            }
            /******************************
             *   MsgQ no-op and comp no-op
             ******************************/
            else if ((ts->input.MsgQ_noop_flag == 1) && (ts->input.comp_noop_flag == 1)){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z;
                  while (row_counts_j > 0){
                     int get_flag = 1;
                     if (get_flag == 1){
                        row_counts_j--;
                     }
                  }
                  int low = T.start[j];
                  int high = T.start[j+1];
                  for (int kk = low; kk < high; kk++){
                     int i = T.i[kk];
                  }
                  j_loc++;
               }
            }
           /*****************
            *   MsgQ no-op
            *****************/
            else if (ts->input.MsgQ_noop_flag == 1){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z;
                  while (row_counts_j > 0){
                     int get_flag = 1;
                     if (get_flag == 1){
                        x_loc[j_loc] -= z;
                        row_counts_j--;
                     }
                  }
                  /* for row j, update elements of x and row_counts */
                  x_loc[j_loc] /= T.diag[j];
                  double xj = x_loc[j_loc];
                  int low = T.start[j];
                  int high = T.start[j+1];
                  for (int kk = low; kk < high; kk++){
                     int i = T.i[kk];
                     z = T.data[kk] * xj;
                  }
                  j_loc++;
               }
            }
           /*****************
            *   comp no-op
            *****************/
            else if (ts->input.comp_noop_flag == 1){
               #pragma omp for schedule(static, lump) nowait
               for (int j = 0; j < n; j++){ /* loop over rows */
                  int row_counts_j = row_counts[j];
                  double z = 0.0;
                  /* stay idle until element j of x is ready to be used */
                  while (row_counts_j > 0){
                     int get_flag = qGet(&Q, j, &z);
                     if (get_flag == 1){
                        row_counts_j--;
                     }
                  }
                  int low = T.start[j];
                  int high = T.start[j+1];
                  for (int kk = low; kk < high; kk++){
                     int i = T.i[kk];
                     qPut(&Q, i, z);
                  }
                  j_loc++;
               }
            }
           /**********************************************
            * standard scheme (no timers, no-ops, etc...) 
            **********************************************/
            else {
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
         }





         /**********************
          *
          *      MsgQ CSR
          *
          **********************/
         else {

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
                  for (int j = 0; j < n; j++){ /* loop over rows */
                     k = my_nzs[kk];
                     int k_end = T.start[j+1];
                     if (k < k_end){ /* does this thread own non-zeros for this row? */
                        int row_counts_j;
                        /* stay idle until element j of x is ready to be used */
                        do {
                           #pragma omp atomic read
                           row_counts_j = row_counts[j];
                        } while (row_counts_j > 0);
                        double xj = x[j];
                        /* for row j, update elements of x and row_counts */
                        while (k < k_end){
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
                  for (int j = 0; j < n; j++){ /* loop over rows */
                     k = my_nzs[kk];
                     int k_end = T.start[j+1];
                     if (k < k_end){ /* does this thread own non-zeros for this row? */
                        int row_counts_j;
                        /* stay idle until element j of x is ready to be used */
                        do {
                           //#pragma omp atomic read
                           row_counts_j = row_counts[j];
                        } while (row_counts_j > 0);
                        double xj = x[j];
                        /* for row j, update elements of x and row_counts */
                        while (k < k_end){
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
                  #pragma omp for schedule(static, lump) nowait
                  for (int i = 0; i < n; i++){ /* loop over rows */
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
                     i_loc++;
                  }
               }
               /************
                * no atomic
                ************/
               else {
                  #pragma omp for schedule(static, lump) nowait
                  for (int i = 0; i < n; i++){ /* loop over rows */
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
                     i_loc++;
                  }
               }
            }
         }
      }


      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;

      if (ts->input.MsgQ_flag == 1){
         i_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            x[i] = x_loc[i_loc];
            i_loc++;
         }
         free(x_loc);

         if (ts->input.MsgQ_wtime_flag == 1){
            ts->output.MsgQ_wtime_vec[tid] = MsgQ_wtime;
         }
         else if (ts->input.MsgQ_cycles_flag == 1){
            ts->output.MsgQ_cycles_vec[tid] = MsgQ_cycles;
         }
         else if (ts->input.comp_wtime_flag == 1){
            ts->output.comp_wtime_vec[tid] = comp_wtime;
         }
         else if (ts->input.comp_cycles_flag == 1){
            ts->output.comp_cycles_vec[tid] = comp_cycles;
         }
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
            }
         }
      }
   }

   if (ts->input.MsgQ_flag == 1){
      //qPrintDAG(&Q);
      qDestroyLock(&Q);
      qFree(&Q);
   }
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
