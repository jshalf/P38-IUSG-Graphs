#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

void TriSolve_FineGrained_COO(TriSolveData *ts,
                              CSR L,
                              CSR U,
                              double *x,
                              double *y,
                              double *b)
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

   int *L_row_counts = (int *)calloc(L_n, sizeof(int));
   int *L_prev_row_counts = (int *)calloc(L_n, sizeof(int));
   int *U_row_counts = (int *)calloc(U_n, sizeof(int));
   int *U_prev_row_counts = (int *)calloc(U_n, sizeof(int));

   int *L_done_flags = (int *)calloc(L_nnz, sizeof(int));
   int *U_done_flags = (int *)calloc(U_nnz, sizeof(int));
   int *U_init_flags = (int *)calloc(U_n, sizeof(int));

   for (int i = 0; i < L_n; i++) L_row_counts[i] = L.i_ptr[i+1] - L.i_ptr[i];
   for (int i = 0; i < U_n; i++) U_row_counts[i] = U.i_ptr[i+1] - U.i_ptr[i] + 1;

   int lump = 1;

   int q_size;
   Queue Q;
   if (ts->input.msgQ_flag == 1){
      Q.type = Q_STDQUEUE;
      q_size = 2*L_n + 2*U_n + L_n + U_n;
      if (Q.type == Q_ARRAY){ 
         int *q_lens = (int *)malloc(q_size * sizeof(int));
         int kk = 0;
         for (int i = 0; i < L_n; i++, kk++) q_lens[kk] = 2;
         for (int i = 0; i < L_n; i++, kk++) q_lens[kk] = 2;
         for (int i = 0; i < U_n; i++, kk++) q_lens[kk] = 2;
         for (int i = 0; i < U_n; i++, kk++) q_lens[kk] = 2;
         for (int i = 0; i < L_n; i++, kk++) q_lens[kk] = L_row_counts[i]+1;
         for (int i = 0; i < U_n; i++, kk++) q_lens[kk] = U_row_counts[i]+1;
         qAlloc(&Q, q_size, q_lens);
         free(q_lens);
      }
      else {
         qAlloc(&Q, q_size, NULL);
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

      #pragma omp for schedule(static, lump)
      for (int i = 0; i < L_n; i++){
         x[i] = b[i] / L.diag[i];
         y[i] = 0.0;
         U_init_count++;
         L_n_loc++;
      }
      #pragma omp for schedule(static, lump)
      for (int i = 0; i < U_n; i++){
         U_n_loc++;
      }
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < L_nnz; k++){
         L_count++;
      }
      L_nnz_loc = L_count;
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < U_nnz; k++){
         U_count++;
      }
      U_nnz_loc = U_count;
      all_count = L_count + U_init_count + U_count;

      CSR L_loc;
      CSR U_loc;
      double *U_diag_loc;
      int *L_done_flags_loc;
      int *U_done_flags_loc;
      int *U_init_flags_loc;
      if ((ts->input.omp_for_flag == 0) || (ts->input.msgQ_flag == 1)){
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
            kk++;
         }
         kk = 0;
         #pragma omp for nowait schedule(static, lump)
         for (int k = 0; k < U_nnz; k++){
            U_loc.i[kk] = U.i[k];
            U_loc.j[kk] = U.j[k];
            U_loc.data[kk] = U.data[k];
            U_loc.diag[kk] = U.diag[U.i[k]];
            kk++;
         }
      }

      int put_stride;
      if (ts->input.msgQ_flag == 1){
         put_stride = 2*L_n + 2*U_n;
              
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < L_n; i++){
            qPut(&Q, i, (double)L_row_counts[i]);
            qPut(&Q, L_n+i, x[i]);
         }
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < U_n; i++){
            int stride = 2*L_n;
            qPut(&Q, stride+i, (double)U_row_counts[i]);
            qPut(&Q, stride+U_n+i, y[i]);
         }
 
         //for (int i = 0; i < Q.size; i++){
         //   printf("%d %f %d %d\n", (i+1) % 10, Q.data[i][0], Q.position[i], Q.max_position[i]);
         //}
         //queue<double> *q = qGetObj();
         //for (int i = 0; i < q_size; i++){
         //   printf("%d\n", q[i].size());
         //}
      }
      #pragma omp barrier
 
      solve_start = omp_get_wtime();
      while (1){ 
         if (ts->input.msgQ_flag == 1){ /* msg Q */
            if (L_count > 0){
               for (int k = 0; k < L_nnz_loc; k++){
                  if (L_done_flags_loc[k] == 0){
                     int j = L_loc.j[k];
                     double L_row_counts_j = 1.0;
                     qPoll(&Q, j, &L_row_counts_j);
                     if ((int)L_row_counts_j == 0){
                        int i = L_loc.i[k];
                        double Lij = L_loc.data[k];
                        double Lii = L_loc.diag[k];

                        double xj, z;
                        qPoll(&Q, L_n + j, &xj);
                        z = -xj * (Lij / Lii);
                        qPut(&Q, put_stride + i, z);

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
                  double L_row_counts_i = 1.0;
                  qPoll(&Q, i, &L_row_counts_i);
                  if (U_init_flags_loc[ii] == 0 && (int)L_row_counts_i == 0){
                     double Uii = U_diag_loc[ii];

                     double xi, z;
                     qPoll(&Q, L_n + i, &xi);
                     z = xi / Uii;
                     qPut(&Q, put_stride+L_n + i, z);

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
                     double U_row_counts_j = 1.0;
                     qPoll(&Q, 2*L_n + j, &U_row_counts_j);
                     if ((int)U_row_counts_j == 0){
                        int i = U_loc.i[k];
                        double Uij = U_loc.data[k];
                        double Uii = U_loc.diag[k];

                        double yj, z;
                        qPoll(&Q, 2*L_n+U_n + j, &yj);
                        z = -yj * (Uij / Uii);
                        qPut(&Q, put_stride+L_n + i, z);

                        U_count--;
                        all_count--;

                        U_done_flags_loc[k] = 1;
                     }
                  }
               }
            }
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < L_n; i++){
               double z;
               while(qGet(&Q, put_stride + i, &z)){
                  qAccum(&Q, L_n+i, z);
                  qAccum(&Q, i, -1.0);
               }
            }
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < U_n; i++){
               double z;
               while(qGet(&Q, put_stride+L_n + i, &z)){
                  qAccum(&Q, 2*L_n+U_n + i, z);
                  qAccum(&Q, 2*L_n + i, -1.0);
               }
            }
         }
         else { /* msg Q */
            if (ts->input.atomic_flag == 1){ /* atomic */
               if (ts->input.omp_for_flag == 1){ /*omp for */ // TODO: fix deadlock when using OpenMP schedule this is not static
                  if (L_count > 0){
                     #pragma omp for schedule(TRISOLVE_OMPFOR_SCHED, lump) nowait
                     for (int k = 0; k < L_nnz; k++){
                        if (L_done_flags[k] == 0){
                           int j = L.j[k];

                           int L_row_counts_j;
                           #pragma omp atomic read
                           L_row_counts_j = L_row_counts[j];
                           if (L_row_counts_j == 0){
                              int i = L.i[k];
                              double Lij = L.data[k];
                              double Lii = L.diag[i];
                              double xj;
                              
                              #pragma omp atomic read
                              xj = x[j];
                              #pragma omp atomic
                              x[i] -= xj * (Lij / Lii);
                              #pragma omp atomic
                              L_row_counts[i]--;

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
                        #pragma omp atomic read
                        L_row_counts_i = L_row_counts[i];
                        if (U_init_flags[i] == 0 && L_row_counts_i == 0){
                           double Uii = U.diag[i];
                           double xi;

                           #pragma omp atomic read
                           xi = x[i];
                           #pragma omp atomic
                           y[i] += xi / Uii;
                           #pragma omp atomic
                           U_row_counts[i]--;

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
               } /* omp for */
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
            } /* atomic */
            else {
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
               } /* omp for */
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
               } /* omp for */
            } /* atomic */
         }/* msg Q */

         ts->output.num_relax[tid]++;

         if (ts->input.async_flag == 0){
            if (all_count == 0){
               converge_flags[tid] = 1;
            }
            #pragma omp barrier
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
            if (all_converge_flag == 1){
               break;
            }
         }
         else {
            if (all_count == 0){
               break;
            }
         }
      }
      ts->output.solve_wtime_vec[tid] = omp_get_wtime() - solve_start;
      #pragma omp barrier

      if (ts->input.msgQ_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < L_n; i++){
            double z;
            qGet(&Q, L_n + i, &z);
            x[i] = z;
         }
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < U_n; i++){
            double z;
            qGet(&Q, 2*L_n+U_n + i, &z);
            y[i] = z;
         }
      }
   }

   free(L_row_counts);
   free(L_prev_row_counts);
   free(U_row_counts);
   free(U_prev_row_counts);

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
}
