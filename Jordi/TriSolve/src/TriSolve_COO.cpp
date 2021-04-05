#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/********************************************************************
 * Async TriSolve (sparse matrix must be in coordinate (COO) format)
 ********************************************************************/
void TriSolve_Async_COO(TriSolveData *ts,
                        CSR T, /* triangular matrix */
                        double *x, /* solution (output) */
                        double *b /* right-hand side */
                        )
{
   int nnz = T.nnz;
   int n = T.n;

   int *row_counts;

   int lump = 1;

   /* put targets are arrays of destination queue ids */ 
   vector<vector<int>> put_targets(n, vector<int>(0));
   int q_size;
   Queue Q;
   /* setup message queue data if message queues are to be used */
   if (ts->input.MsgQ_flag == 1){
      q_size = nnz + n;
      qAlloc(&Q, q_size);

      int get_stride = 0;
      for (int k = 0; k < nnz; k++){
         int j = T.j[k];
         put_targets[j].push_back(get_stride+k);
      }
      get_stride += nnz;
      for (int i = 0; i < n; i++){
         put_targets[i].push_back(get_stride+i);
      }
      qInitLock(&Q);
   }
   else {
      /* row counts to track updates to solution vectors */
      row_counts = (int *)calloc(n, sizeof(int));
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = ts->input.num_threads;
      int accum_count = 0;
      double solve_start, setup_start;
      int nnz_loc = 0;
      int n_loc = 0;
      int kk, ii, m;
      CSR T_loc;
      int *my_nnzs, *my_rows;
      int *done_flags_loc;
      int *row_counts_loc;

      setup_start = omp_get_wtime();

      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n; i++){
         n_loc++;
      }
      #pragma omp for schedule(static, lump) nowait
      for (int k = 0; k < nnz; k++){
         nnz_loc++;
      }
      /* number of non-zeros to loop over */
      int block_size;// = std::max(1, std::min(ts->input.block_size, nnz_loc));
      /* accum_count is the total number of updates done by this thread.  
       * when accum_count is decrimented to zero, convergence is achieved. */ 
      accum_count = nnz_loc;

      //T_loc.j = (int *)calloc(nnz_loc, sizeof(int));
      //T_loc.i = (int *)calloc(nnz_loc, sizeof(int));
      //T_loc.data = (double *)calloc(nnz_loc, sizeof(double));
      //T_loc.diag = (double *)calloc(nnz_loc, sizeof(double));
      my_nnzs = (int *)calloc(nnz_loc, sizeof(int));
      done_flags_loc = (int *)calloc(nnz_loc, sizeof(int));
      kk = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int k = 0; k < nnz; k++){
         //T_loc.j[kk] = T.j[k];
         //T_loc.i[kk] = T.i[k];
         //T_loc.data[kk] = T.data[k];
         //T_loc.diag[kk] = T.diag[T.i[k]];
         my_nnzs[kk] = k;
         kk++;
      }

      //int num_blocks = nnz_loc / block_size;
      //int *kk_blocks = (int *)calloc(num_blocks+1, sizeof(int));
      //for (int m = 0; m < num_blocks; m++){
      //   if (m == num_blocks-1){
      //      kk_blocks[m+1] = nnz_loc;
      //   }
      //   else {
      //      kk_blocks[m+1] = kk_blocks[m] + block_size;
      //   }
      //}

      /* compute row counts (number of non-zeros per row) and other data */
      if (ts->input.MsgQ_flag == 1){
         accum_count += n_loc;
         my_rows = (int *)calloc(n_loc, sizeof(int));
         row_counts_loc = (int *)calloc(n_loc, sizeof(int));
         ii = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            my_rows[ii] = i;
            row_counts_loc[ii] = T.i_ptr[i+1] - T.i_ptr[i];
            ii++;
         }
      }
      else {
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n; i++){
            row_counts[i] = T.i_ptr[i+1] - T.i_ptr[i];
         }
      }
      /* these strides denote starting points in the array of queues for a thread to get from or put to */
      int get_stride, put_stride;

      ts->output.setup_wtime_vec[tid] = omp_get_wtime() - setup_start;
      solve_start = omp_get_wtime();

      /* initialize solution */
      #pragma omp for schedule(static, lump)
      for (int i = 0; i < n; i++){
         x[i] = b[i] / T.diag[i];
      }

      //solve_start = omp_get_wtime();
      int num_relax = 0, num_iters = 0;
      /* iterate until convergence is achieved */
      /******************************* 
       * MESSAGE QUEUE IMPLEMENTATION
       *******************************/
      if (ts->input.MsgQ_flag == 1){
         int kk_start = 0;
         int kk_end = block_size;
         int ii_start = 0;
         int ii_end = block_size;
         int kk_block_count = block_size;
         int ii_block_count = block_size;
         while (1){
            get_stride = 0;
            put_stride = nnz;
            for (int kk = kk_start; kk < kk_end; kk++){ /* loop over block of non-zeros */
               int k;
               if (done_flags_loc[kk] == 0){ /* is the update for this non-zero done? */
                  k = my_nnzs[kk];
                  double xj;
                  if (qGet(&Q, get_stride+k, &xj)){ /* check if element j of x has arrived */
                     int i = T.i[k];
                     int j = T.j[k];
                     double Tij = T.data[k];
                     double Tii = T.diag[i];

                     /* compute the update to x */
                     double z = -xj * (Tij / Tii);
                     /* send update using a put primitive */
                     qPut(&Q, put_stride+i, z);

                     /* decriment counters */
                     accum_count--;
                     if (accum_count == 0){
                        break;
                     }

                     /* update block of non-zeros */
                     kk_end = std::min(kk_end+1, nnz_loc);
                     if (kk == kk_start){
                        kk_start++;
                     }

                     //kk_block_count--;

                     done_flags_loc[kk] = 1;

                     num_relax++;
                  }
               }
            }
            if (accum_count == 0){
               break;
            }
           // if (kk_block_count == 0){
           //    kk_block_count = block_size;
           //    kk_start = kk_end;
           //    kk_end = std::min(kk_end+block_size, nnz_loc);
           // }

            /* receive updates to x */
            get_stride = nnz;
            for (int ii = ii_start; ii < ii_end; ii++){ /* loop over block of rows */
               int i;
               if (row_counts_loc[ii] > 0){ /* if there are no more messages to get */
                  i = my_rows[ii];
                  double z;
                  while (qGet(&Q, get_stride+i, &z)){ /* get update */
                     x[i] += z; /* increment x*/
                     row_counts_loc[ii]--; /* decrement row counts */
                  }
               }
               /* send x[i] to rows that need it */
               if (row_counts_loc[ii] == 0){ /* if this update is the last */
                  i = my_rows[ii];
                  double xi = x[i];
                  for (int p = 0; p < put_targets[i].size(); p++){
                     qPut(&Q, put_targets[i][p], xi); /* send completed element of x to queues that need it */
                  }
                  accum_count--;
                  if (accum_count == 0){
                     break;
                  }
                  //ii_block_count--;
                  /* update block  of rows */
                  ii_end = std::min(ii_end+1, n_loc);
                  if (ii == ii_start){
                     ii_start++;
                  }
                  row_counts_loc[ii]--;
               }
            }
            if (accum_count == 0){
               break;
            }
           // if (ii_block_count == 0){
           //    ii_block_count = block_size;
           //    ii_start = ii_end;
           //    ii_end = std::min(ii_end+block_size, n_loc);
           // }

            //printf("%d %d %d %d %d %d\n", ii_start, kk_start, ii_end, kk_end, ii_block_count, kk_block_count);
         }
      }
      else {
         if (ts->input.atomic_flag == 1){ /* atomic implementation using atomics */
            /*******************************************************************
             * ATOMICS IMPLEMENTATION 
             *******************************************************************/
            int block_count = block_size;
            int kk_start = 0, kk_end = block_size;
            m = 0;
            while (1){
               //for (int kk = kk_blocks[m]; kk < kk_blocks[m+1]; kk++){
               for (int kk = kk_start; kk < kk_end; kk++){ /* loop over block of non-zeros */
                  int k;
                  if (done_flags_loc[kk] == 0){ /* is the update for this non-zero done? */
                     k = my_nnzs[kk];
                     int j = T.j[k];
                     int row_counts_j;

                     /* check if x[j] is available */
                     #pragma omp atomic read
                     row_counts_j = row_counts[j];

                     if (row_counts_j == 0){
                        int i = T.i[k];
                        double Tij = T.data[k];
                        double Tii = T.diag[i];

                        /* update x and row counts for i*/
                        #pragma omp atomic
                        x[i] -= x[j] * (Tij / Tii);
                        #pragma omp atomic
                        row_counts[i]--;

                        accum_count--;
                        if (accum_count == 0){
                           break;
                        }

                        //block_count--;
                        /* update block of non-zeros */
                        kk_end = std::min(kk_end+1, nnz_loc);
                        if (kk == kk_start){
                           kk_start++;
                        }

                        done_flags_loc[kk] = 1;

                        num_relax++;
                     }
                  }
               }
               if (accum_count == 0){
                  break;
               }
               //if (block_count == 0){
               //   block_count = block_size;
               //   m++;
               //}
            }
         }
         else {
            /**************************************************************
             * NO-ATOMICS IMPLEMENTATION
             * (this implementation is the same as the atomics implementation 
             *  but with atomics removed. the atomics are removed for 
             *  performance analysis and therefore
             *  this implementation will not produce a correct result.
             **************************************************************/
            int block_count = block_size;
            int kk_start = 0, kk_end = block_size;
            m = 0;
            while (1){ 
               //for (int kk = kk_blocks[m]; kk < kk_blocks[m+1]; kk++){
               for (int kk = kk_start; kk < kk_end; kk++){ /* loop over block of non-zeros */
                  int k;
                  if (done_flags_loc[kk] == 0){ /* is the update for this non-zero done? */
                     k = my_nnzs[kk];
                     int j = T.j[k];
                     int row_counts_j;

                     /* check if x[j] is available */
                     //#pragma omp atomic read
                     row_counts_j = row_counts[j];

                     if (row_counts_j == 0){
                        int i = T.i[k];
                        double Tij = T.data[k];
                        double Tii = T.diag[i];

                        /* update x and row counts for i*/
                        //#pragma omp atomic
                        x[i] -= x[j] * (Tij / Tii);
                        #pragma omp atomic
                        row_counts[i]--;

                        accum_count--;
                        if (accum_count == 0){
                           break;
                        }

                        //block_count--;
                        /* update block of non-zeros */
                        kk_end = std::min(kk_end+1, nnz_loc);
                        if (kk == kk_start){
                           kk_start++;
                        }

                        done_flags_loc[kk] = 1;

                        num_relax++;
                     }
                  }
               }
               if (accum_count == 0){
                  break;
               }
               //if (block_count == 0){
               //   block_count = block_size;
               //   m++;
               //}
            }
         }
      }

      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;

      ts->output.num_relax[tid] = num_relax;
      ts->output.num_iters[tid] = num_iters;

      //free(T_loc.j);
      //free(T_loc.i);
      //free(T_loc.diag);
      //free(T_loc.data);
      free(my_nnzs);
      free(done_flags_loc);
      //free(kk_blocks);
      if (ts->input.MsgQ_flag == 1){
         free(my_rows);
         free(row_counts_loc);
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
