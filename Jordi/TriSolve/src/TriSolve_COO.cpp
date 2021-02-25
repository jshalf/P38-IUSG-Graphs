#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/*************************************************************************************************
 * ITERATIVE FINE-GRAINED ASYNCHRONOUS TRISOLVE:
 *    Sparse matrix must be in coordinate (COO) format.
 *    Each thread is assigned a partion of the rows and non-zeros of L and U (four partiions).
 *    Thread t iterates until all non-zeros of L and U in its partition have been used.
 *    During each iteration, thread t carries out three general steps:
 *     1. Loop over t's non-zeros of L and compute xi -= xj Lij / Lii if xj is available
 *     2. Loop over t's rows of U and compute yi = xi / Uii if xi is available
 *     3. Loop over t's non-zeros of U and compute yi -= yj Uij / Uii if yj is available
 *    Initially x = b ./ diag(L) before iterating.
 *************************************************************************************************/

void TriSolve_FineGrained_COO(TriSolveData *ts,
                              CSR L, /* lower triangular part of matrix */
                              CSR U, /* upper triangular part of matrix */
                              double *x, /* solution for lower tri-solve (output) */
                              double *y, /* solution for upper tri-solve (output) */
                              double *b /* right-hand side */
                              )
{
   int L_nnz = L.nnz;
   int U_nnz = U.nnz;
   int L_n = L.n;
   int U_n = U.n;

   int *converge_flags = (int *)calloc(ts->input.num_threads, sizeof(int));
   int all_converge_flag = 0;

   int *dummy_L_row_counts, *dummy_U_row_counts;
   double *dummy_x, *dummy_y;
   
   if (ts->input.atomic_flag == 0){
      dummy_L_row_counts = (int *)calloc(L_n, sizeof(int));
      dummy_U_row_counts = (int *)calloc(U_n, sizeof(int));
      dummy_x = (double *)calloc(L_n, sizeof(double));
      dummy_y = (double *)calloc(L_n, sizeof(double));
   }

   /* row counts to track updates to solution vectors */
   int *L_row_counts = (int *)calloc(L_n, sizeof(int));
   int *U_row_counts = (int *)calloc(U_n, sizeof(int));

   /* flags to indicate when a trisolve is complete */
   int *L_done_flags = (int *)calloc(L_nnz, sizeof(int)); /* lower */
   int *U_done_flags = (int *)calloc(U_nnz, sizeof(int)); /* upper */
   int *U_init_flags = (int *)calloc(U_n, sizeof(int)); /* initialization of upper */

   /* initial row counts are the number of non-zeros in each row.  upon convergence, row counts will be zero */
   for (int i = 0; i < L_n; i++) {
      L_row_counts[i] = L.i_ptr[i+1] - L.i_ptr[i];
   }
   for (int i = 0; i < U_n; i++) {
      U_row_counts[i] = U.i_ptr[i+1] - U.i_ptr[i] + 1;
   }

   int lump = 1;

   /* put targets are arrays of destination queue ids */ 
   vector<vector<int>> L_put_targets(L_n, vector<int>(0)), U_put_targets(U_n, vector<int>(0));
   int q_size;
   Queue Q;
   /* setup message queue data if message queues are to be used */
   if (ts->input.MsgQ_flag == 1){
      q_size = L_nnz + U_n + L_nnz + L_n + U_n;
      qAlloc(&Q, q_size);

      int get_stride = 0;
      for (int k = 0; k < L_nnz; k++){
         int j = L.j[k];
         L_put_targets[j].push_back(get_stride+k);
      }
      get_stride += L_nnz;
      for (int i = 0; i < U_n; i++){
         L_put_targets[i].push_back(get_stride+i);
      }
      get_stride += U_n;
      for (int k = 0; k < U_nnz; k++){
         int j = U.j[k];
         U_put_targets[j].push_back(get_stride+k);
      }

      qInitLock(&Q);
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = ts->input.num_threads;
      int all_count = 0, L_count = 0, U_init_count = 0, U_count = 0;
      double solve_start, start;
      int L_nnz_loc = 0;
      int U_nnz_loc = 0;
      int U_n_loc = 0;
      int L_n_loc = 0;
      int kk, ii;

      /* initialize solutions */
      #pragma omp for schedule(static, lump)
      for (int i = 0; i < L_n; i++){
         x[i] = b[i] / L.diag[i];
         y[i] = 0.0;
         L_n_loc++; /* number of L rows for this thread */
      }
      #pragma omp for schedule(static, lump)
      for (int i = 0; i < U_n; i++){
         U_n_loc++; /* number of U rows for this thread */
      }
      U_init_count = U_n_loc;
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < L_nnz; k++){
         L_count++; /* number of L nnz for this thread */
      }
      L_nnz_loc = L_count;
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < U_nnz; k++){
         U_count++; /* number of U nnz for this thread */
      }
      U_nnz_loc = U_count;
      /* all_count is the total number of updates done by this thread.  all_count is decrimented to zero, convergence is achieved. */ 
      all_count = L_count + U_init_count + U_count;

      CSR L_loc;
      CSR U_loc;
      double *U_diag_loc;
      int *L_done_flags_loc, *U_done_flags_loc;
      int *U_init_flags_loc;
      int *L_nnz_map, *U_nnz_map;
      /* if we are not using OpenMP for loops or we are using message queues, make some data local. */ 
      if ((ts->input.omp_for_flag == 0) || (ts->input.MsgQ_flag == 1)){
         U_loc.i = (int *)calloc(U_nnz_loc, sizeof(int));
         U_loc.j = (int *)calloc(U_nnz_loc, sizeof(int));
         U_loc.data = (double *)calloc(U_nnz_loc, sizeof(double));
         U_loc.diag = (double *)calloc(U_nnz_loc, sizeof(double));
         L_loc.i = (int *)calloc(L_nnz_loc, sizeof(int));
         L_loc.j = (int *)calloc(L_nnz_loc, sizeof(int));
         L_loc.data = (double *)calloc(L_nnz_loc, sizeof(double));
         L_loc.diag = (double *)calloc(L_nnz_loc, sizeof(double));
         U_diag_loc = (double *)calloc(U_n_loc, sizeof(double));

         L_done_flags_loc = (int *)calloc(L_nnz_loc, sizeof(int));
         U_done_flags_loc = (int *)calloc(U_nnz_loc, sizeof(int));
         U_init_flags_loc = (int *)calloc(U_n_loc, sizeof(int));

         L_nnz_map = (int *)calloc(L_nnz_loc, sizeof(int));
         U_nnz_map = (int *)calloc(U_nnz_loc, sizeof(int));         

         ii = 0;
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < U_n; i++){
            U_diag_loc[ii] = U.diag[i];
            ii++;
         }

         kk = 0;
         #pragma omp for nowait schedule(static, lump)
         for (int k = 0; k < L_nnz; k++){
            L_loc.i[kk] = L.i[k];
            L_loc.j[kk] = L.j[k];
            L_loc.data[kk] = L.data[k];
            L_loc.diag[kk] = L.diag[L.i[k]];
            L_nnz_map[kk] = k;
            kk++;
         }
         kk = 0;
         #pragma omp for nowait schedule(static, lump)
         for (int k = 0; k < U_nnz; k++){
            U_loc.i[kk] = U.i[k];
            U_loc.j[kk] = U.j[k];
            U_loc.data[kk] = U.data[k];
            U_loc.diag[kk] = U.diag[U.i[k]];
            U_nnz_map[kk] = k;
            kk++;
         }
      }
      if (ts->input.MsgQ_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < L_n; i++){
            all_count++;
         }
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < U_n; i++){
            all_count++;
         }
      }
      /* these strides denote starting points in the array of queues for a thread to get from or put to */
      int get_stride, put_stride;
      #pragma omp barrier
 
      /* iterate until convergence is achieved */
      solve_start = omp_get_wtime();
      while (1){
         /************************************************************** 
          * MESSAGE QUEUE IMPLEMENTATION
          **************************************************************/
         if (ts->input.MsgQ_flag == 1){
            /* intialize strides */
            get_stride = 0;
            put_stride = L_nnz + U_n + U_nnz;
            /* compute Lx=b */
            if (L_count > 0){ /* if solution of Lx=b is not complete */
               for (int k = 0; k < L_nnz_loc; k++){ /* loop over non-zero values of L */
                  if (L_done_flags_loc[k] == 0){ /* if this non-zero has not yet been used to update x */
                     double xj;
                     if (qGet(&Q, get_stride+L_nnz_map[k], &xj)){ /* check if element j or x has arrived */
                        int i = L_loc.i[k];
                        int j = L_loc.j[k];
                        double Lij = L_loc.data[k];
                        double Lii = L_loc.diag[k];
                        
                        /* compute the update to x */
                        double z = -xj * (Lij / Lii);
                        /* send update using a put primitive */
                        qPut(&Q, put_stride+i, z);

                        /* decriment counters */
                        L_count--;
                        all_count--;

                        /* the update from this non-zero is now complete */
                        L_done_flags_loc[k] = 1;
                     }
                  }
               }
            }
            /* update strides */
            get_stride += L_nnz;
            put_stride += L_n;
            /* next two if statements are for computing Uy=x */
            if (U_init_count > 0){ /* if initialization of solution Uy=x is not complete (initialization requires solution from Lx=b) */
               ii = 0;
               #pragma omp for schedule(static, lump) nowait /* loop over rows of U */
               for (int i = 0; i < U_n; i++){
                  if (U_init_flags_loc[ii] == 0){ /* if this row has has not yet been initialized */ 
                     double xi;
                     if (qGet(&Q, get_stride+i, &xi) == 1){ /* check if element i or x has arrived */
                        double Uii = U_diag_loc[ii];

                        /* compute the update to y */
                        double z = xi / Uii;
                        /* send update using a put primitive */
                        qPut(&Q, put_stride+i, z);

                        /* decriment counters */
                        U_init_count--;
                        all_count--;

                        /* the initialization of this row is now complete */
                        U_init_flags_loc[ii] = 1;
                     }
                  }
                  ii++;
               }
            }
            get_stride += U_n;
            if (U_count > 0){ /* if solution of Uy=x is not complete (see Lx=b block above for similar comments) */
               for (int k = 0; k < U_nnz_loc; k++){
                  if (U_done_flags_loc[k] == 0){
                     double yj;
                     if (qGet(&Q, get_stride+U_nnz_map[k], &yj)){
                        int i = U_loc.i[k];
                        int j = U_loc.j[k];
                        double Uij = U_loc.data[k];
                        double Uii = U_loc.diag[k];

                        double z = -yj * (Uij / Uii);
                        qPut(&Q, put_stride+i, z);

                        U_count--;
                        all_count--;

                        U_done_flags_loc[k] = 1;
                     }
                  }
               }
            }

            /* receive updates to x */
            get_stride = L_nnz + U_n + U_nnz;
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < L_n; i++){
               if (L_row_counts[i] > 0){ /* if there are no more messages to get */
                  double z;
                  while (qGet(&Q, get_stride+i, &z)){ /* get update */
                     x[i] += z; /* increment x*/
                     L_row_counts[i]--; /* decrement row counts */
                  }
               }
               if (L_row_counts[i] == 0){ /* if this update is the last */
                  double xi = x[i];
                  for (int p = 0; p < L_put_targets[i].size(); p++){
                     qPut(&Q, L_put_targets[i][p], xi); /* send completed element of x to queues that need it */
                  }
                  /* decrement counters */
                  L_row_counts[i]--;
                  all_count--;
               }
            }
            /* receive updates to y (for details, see above comments regarding x) */
            get_stride += L_n;
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < U_n; i++){
               if (U_row_counts[i] > 0){
                  double z;
                  while (qGet(&Q, get_stride+i, &z)){
                     y[i] += z;
                     U_row_counts[i]--;
                  }
               }
               if (U_row_counts[i] == 0){
                  double yi = y[i];
                  for (int p = 0; p < U_put_targets[i].size(); p++){
                     qPut(&Q, U_put_targets[i][p], yi);
                  }
                  U_row_counts[i]--;
                  all_count--;
               }
            }
         }
         else {
            if (ts->input.atomic_flag == 1){ /* atomic implementation using atomics */
               /**************************************************************
                * ATOMICS IMPLEMENTATION USING OPENMP FOR LOOPS
                **************************************************************/
               if (ts->input.omp_for_flag == 1){ /* atomic impelementation using omp for */ // TODO: fix deadlock when using OpenMP non-static schedule
                  /* compute Lx=b */
                  if (L_count > 0){ /* if solution of Lx=b is not complete */
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait /* loop over non-zero values of L */
                     for (int k = 0; k < L_nnz; k++){
                        if (L_done_flags[k] == 0){ /* if this non-zero has not yet been used to update x */
                           int j = L.j[k];
                           int L_row_counts_j;
                           #pragma omp atomic read /* atomically read update counter for element j */
                           L_row_counts_j = L_row_counts[j]; 
                           if (L_row_counts_j == 0){ /* if L_row_counts[j] is zero, then all updates for element j are complete */
                              int i = L.i[k];
                              double Lij = L.data[k];
                              double Lii = L.diag[i];
                              double xj;
                              
                              /* do the update for this non-zero using atomics */
                              #pragma omp atomic read /* atomically read element j */
                              xj = x[j];
                              #pragma omp atomic /* atomically update element i */
                              x[i] -= xj * (Lij / Lii);
                              #pragma omp atomic /* atomically decrement the counter for element i */
                              L_row_counts[i]--;

                              /* decriment counters */
                              L_count--;
                              all_count--;

                              /* the update from this non-zero is now complete */
                              L_done_flags[k] = 1;
                           }
                        }
                     }
                  }
                  /* next two if statements are for computing Uy=x */
                  if (U_init_count > 0){ /* if initialization of solution Uy=x is not complete (initialization requires solution from Lx=b) */
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait /* loop over rows of U */
                     for (int i = 0; i < U_n; i++){
                        int L_row_counts_i;
                        #pragma omp atomic read /* atomically read update counter for element i of x */
                        L_row_counts_i = L_row_counts[i];
                        if (U_init_flags[i] == 0 && L_row_counts_i == 0){ /* if element i of y has not been initialized AND element i of x is not done updating */
                           double Uii = U.diag[i];
                           double xi;

                           /* compute update to element i of y using atomics */
                           #pragma omp atomic read /* atomically read element i of x */
                           xi = x[i];
                           #pragma omp atomic /* atomically update element i of y */
                           y[i] += xi / Uii;
                           #pragma omp atomic /* atomically decrement the counter for element i */
                           U_row_counts[i]--;

                           /* decriment counters */
                           U_init_count--;
                           all_count--;

                           /* the initialization of element i of y is now complete */
                           U_init_flags[i] = 1;
                        }
                     }
                  }
                  if (U_count > 0){ /* if solution of Uy=x is not complete (see Lx=b block above for similar comments) */
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait
                     for (int k = 0; k < U_nnz; k++){
                        if (U_done_flags[k] == 0){
                           int j = U.j[k];
                           int U_row_counts_j;
                           #pragma omp atomic read
                           U_row_counts_j = U_row_counts[j];
                           if (U_row_counts_j == 0){
                              int i = U.i[k];
                              double Uij = U.data[k];
                              double Uii = U.diag[i];
                              double yj;
                              
                              #pragma omp atomic read
                              yj = y[j];
                              #pragma omp atomic
                              y[i] -= yj * (Uij / Uii);
                              #pragma omp atomic
                              U_row_counts[i]--;

                              U_count--;
                              all_count--;

                              U_done_flags[k] = 1;
                           }
                        }
                     }
                  }
               }
               /*******************************************************************
                * ATOMICS IMPLEMENTATION WITHOUT OPENMP FOR LOOPS 
                * (see ``atomics implementation using OpenMP for loops'' 
                *  implementation block above for detailed comments)
                *******************************************************************/
               else {
                  if (L_count > 0){
                     for (int k = 0; k < L_nnz_loc; k++){
                        if (L_done_flags_loc[k] == 0){
                           int j = L_loc.j[k];
                           int L_row_counts_j;
                           #pragma omp atomic read
                           L_row_counts_j = L_row_counts[j];
                           if (L_row_counts_j == 0){
                              int i = L_loc.i[k];
                              double Lij = L_loc.data[k];
                              double Lii = L_loc.diag[k];
                              double xj;

                              #pragma omp atomic read
                              xj = x[j];
                              #pragma omp atomic
                              x[i] -= xj * (Lij / Lii);
                              #pragma omp atomic
                              L_row_counts[i]--;

                              L_count--;
                              all_count--;

                              L_done_flags_loc[k] = 1;
                           }
                        }
                     }
                  }
                  if (U_init_count > 0){
                     ii = 0;
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < U_n; i++){
                        int L_row_counts_i;
                        #pragma omp atomic read
                        L_row_counts_i = L_row_counts[i];
                        if (U_init_flags_loc[ii] == 0 && L_row_counts_i == 0){
                           double Uii = U_diag_loc[ii];
                           double xi;
                           
                           #pragma omp atomic read
                           xi = x[i];
                           #pragma omp atomic
                           y[i] += xi / Uii;
                           #pragma omp atomic
                           U_row_counts[i]--;

                           U_init_flags_loc[ii] = 1;
                           U_init_count--;
                           all_count--;
                        }
                        ii++;
                     }
                  }
                  if (U_count > 0){
                     for (int k = 0; k < U_nnz_loc; k++){
                        if (U_done_flags_loc[k] == 0){
                           int j = U_loc.j[k];
                           int U_row_counts_j;
                           #pragma omp atomic read
                           U_row_counts_j = U_row_counts[j];
                           if (U_row_counts_j == 0){
                              int i = U_loc.i[k];
                              double Uij = U_loc.data[k];
                              double Uii = U_loc.diag[k];
                              double yj;

                              #pragma omp atomic read
                              yj = y[j];
                              #pragma omp atomic
                              y[i] -= yj * (Uij / Uii);
                              #pragma omp atomic
                              U_row_counts[i]--;

                              U_count--;
                              all_count--;

                              U_done_flags_loc[k] = 1;
                           }
                        }
                     }
                  }
               }
            }
            else {
               /**************************************************************
                * NO-ATOMICS IMPLEMENTATION USING OPENMP FOR LOOPS
                * (this implementation is the same as the ``atomics
                *  implementation using OpenMP for loops'' but with atomics removed.
                *  the atomics are removed for performance analysis.
                *  this implementation will not produce a correct result.
                *  see ``atomics implementation using OpenMP for loops''
                *  implementation block above for detailed comments.)
                **************************************************************/ 
               if (ts->input.omp_for_flag == 1){
                  if (L_count > 0){
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait
                     for (int k = 0; k < L_nnz; k++){
                        if (L_done_flags[k] == 0){
                           int j = L.j[k];
                           int L_row_counts_j;

                           //L_row_counts_j = L_row_counts[j];
                           //start = omp_get_wtime();
                           //#pragma omp atomic read
                           L_row_counts_j = L_row_counts[j];
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                           if (L_row_counts_j == 0){
                              int i = L.i[k];
                              double Lij = L.data[k];
                              double Lii = L.diag[i];

                              //start = omp_get_wtime();
                              //#pragma omp atomic
                              x[i] -= x[j] * (Lij / Lii);
                              #pragma omp atomic
                              L_row_counts[i]--;
                              //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                              //dummy_x[i] -= dummy_x[j] * (Lij / Lii);
                              //dummy_L_row_counts[i]--;

                              L_count--;
                              all_count--;

                              L_done_flags[k] = 1;
                           }
                        }
                     }
                  }
                  if (U_init_count > 0){
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait
                     for (int i = 0; i < U_n; i++){
                        int L_row_counts_i;

                        //L_row_counts_i = L_row_counts[i];
                        //start = omp_get_wtime();
                        //#pragma omp atomic read
                        L_row_counts_i = L_row_counts[i];
                        //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                        if (U_init_flags[i] == 0 && L_row_counts_i == 0){
                           double Uii = U.diag[i];

                           //start = omp_get_wtime();
                           //#pragma omp atomic
                           y[i] += x[i] / Uii;
                           #pragma omp atomic
                           U_row_counts[i]--;
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                           //dummy_y[i] += dummy_x[i] / Uii;
                           //dummy_U_row_counts[i]--;

                           U_init_flags[i] = 1;
                           U_init_count--;
                           all_count--;
                        }
                     }
                  }
                  if (U_count > 0){
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait
                     for (int k = 0; k < U_nnz; k++){
                        if (U_done_flags[k] == 0){
                           int j = U.j[k];
                           int U_row_counts_j;

                           //U_row_counts_j = U_row_counts[j];
                           //start = omp_get_wtime();
                           //#pragma omp atomic read
                           U_row_counts_j = U_row_counts[j];
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                           if (U_row_counts_j == 0){
                              int i = U.i[k];
                              double Uij = U.data[k];
                              double Uii = U.diag[i];

                              //start = omp_get_wtime();
                              //#pragma omp atomic
                              y[i] -= y[j] * (Uij / Uii);
                              #pragma omp atomic
                              U_row_counts[i]--;
                              //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                              //dummy_y[i] -= dummy_y[j] * (Uij / Uii);
                              //dummy_U_row_counts[i]--;

                              U_count--;
                              all_count--;

                              U_done_flags[k] = 1;
                           }
                        }
                     }
                  }
               }
               /**************************************************************
                * NO-ATOMICS IMPLEMENTATION WITHOUT OPENMP FOR LOOPS
                * (this implementation is the same as the ``atomics
                *  implementation without OpenMP for loops'' but with atomics removed.
                *  the atomics are removed for performance analysis.
                *  this implementation will not produce a correct result.
                *  see ``atomics implementation using OpenMP for loops''
                *  implementation block above for detailed comments.)
                **************************************************************/
               else {
                  if (L_count > 0){
                     for (int k = 0; k < L_nnz_loc; k++){
                        if (L_done_flags_loc[k] == 0){
                           int j = L_loc.j[k];
                           int L_row_counts_j;

                           //L_row_counts_j = L_row_counts[j];
                           //start = omp_get_wtime();
                           //#pragma omp atomic read
                           L_row_counts_j = L_row_counts[j];
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                           if (L_row_counts_j == 0){
                              int i = L_loc.i[k];
                              double Lij = L_loc.data[k];
                              double Lii = L_loc.diag[k];

                              //start = omp_get_wtime();
                              //#pragma omp atomic
                              x[i] -= x[j] * (Lij / Lii);
                              #pragma omp atomic
                              L_row_counts[i]--;
                              //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                              //dummy_x[i] -= dummy_x[j] * (Lij / Lii);
                              //dummy_L_row_counts[i]--;

                              L_count--;
                              all_count--;

                              L_done_flags_loc[k] = 1;
                           }
                        }
                     }
                  }
                  if (U_init_count > 0){
                     ii = 0;
                     #pragma omp for schedule(static, lump) nowait
                     for (int i = 0; i < U_n; i++){
                        int L_row_counts_i;

                        //L_row_counts_i = L_row_counts[i];
                        //start = omp_get_wtime();
                        //#pragma omp atomic read
                        L_row_counts_i = L_row_counts[i];
                        //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                        if (U_init_flags_loc[ii] == 0 && L_row_counts_i == 0){
                           double Uii = U_diag_loc[ii];

                           //start = omp_get_wtime();
                           //#pragma omp atomic
                           y[i] += x[i] / Uii;
                           #pragma omp atomic
                           U_row_counts[i]--;
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                           //dummy_y[i] += dummy_x[i] / Uii;
                           //dummy_U_row_counts[i]--;

                           U_init_flags_loc[ii] = 1;
                           U_init_count--;
                           all_count--;
                        }
                        ii++;
                     }
                  }
                  if (U_count > 0){
                     for (int k = 0; k < U_nnz_loc; k++){
                        if (U_done_flags_loc[k] == 0){
                           int j = U_loc.j[k];
                           int U_row_counts_j;

                           //U_row_counts_j = U_row_counts[j];
                           //start = omp_get_wtime();
                           //#pragma omp atomic read
                           U_row_counts_j = U_row_counts[j];
                           //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;

                           if (U_row_counts_j == 0){
                              int i = U_loc.i[k];
                              double Uij = U_loc.data[k];
                              double Uii = U_loc.diag[k];

                              //start = omp_get_wtime();
                              //#pragma omp atomic
                              y[i] -= y[j] * (Uij / Uii);
                              #pragma omp atomic
                              U_row_counts[i]--;
                              //ts->output.atomic_wtime_vec[tid] += omp_get_wtime() - start;
                              //dummy_y[i] -= dummy_y[j] * (Uij / Uii);
                              //dummy_U_row_counts[i]--;

                              U_count--;
                              all_count--;

                              U_done_flags_loc[k] = 1;
                           }
                        }
                     }
                  }
               }
            }
         }

         ts->output.num_relax[tid]++;

         /* check for convergence */
         if (ts->input.async_flag == 0){ /* synchronous version */
            /* if the counter has reached zero, indicate convergence for this thread */
            if (all_count == 0){
               converge_flags[tid] = 1;
            }
            #pragma omp barrier
            /* thread 0 now checks to see if all threads have converged.
             * if true, set a global flag to indicate that convergence has been achieved */
            if (tid == 0){
               int all_converge_flag_loc = 1;
               for (int t = 0; t < ts->input.num_threads; t++){
                  if (converge_flags[t] == 0){
                     all_converge_flag_loc = 0;
                     break;
                  }
               }
               if (all_converge_flag_loc == 1){
                  all_converge_flag = 1;
               }
            }
            #pragma omp barrier
            /* stop iterating if convergence is achieved */
            if (all_converge_flag == 1){
               break;
            }
         }
         else {
            /* in the asynchronous version, a thread stops once the counter reaches 0 */
            if (all_count == 0){
               break;
            }
         }
      }
      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;
      #pragma omp barrier
   }

   free(L_row_counts);
   free(U_row_counts);

   if (ts->input.atomic_flag == 0){
      free(dummy_L_row_counts);
      free(dummy_U_row_counts);
      free(dummy_x);
      free(dummy_y);
   }

   free(L_done_flags);
   free(U_done_flags);
   free(U_init_flags);

   free(converge_flags);

   if (ts->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
