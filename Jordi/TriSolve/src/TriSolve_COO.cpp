#include "TriSolve.hpp"
#include "../../src/Matrix.hpp"

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

   int *dummy_L_row_counts = (int *)calloc(L_n, sizeof(int));
   int *dummy_U_row_counts = (int *)calloc(U_n, sizeof(int));
   double *dummy_x = (double *)calloc(L_n, sizeof(double));
   double *dummy_y = (double *)calloc(L_n, sizeof(double));

   //volatile int *L_row_counts = (int *)calloc(L_n, sizeof(int));
   //volatile int *L_prev_row_counts = (int *)calloc(L_n, sizeof(int));
   //volatile int *U_row_counts = (int *)calloc(U_n, sizeof(int));
   //volatile int *U_prev_row_counts = (int *)calloc(U_n, sizeof(int));
   int *L_row_counts = (int *)calloc(L_n, sizeof(int));
   int *L_prev_row_counts = (int *)calloc(L_n, sizeof(int));
   int *U_row_counts = (int *)calloc(U_n, sizeof(int));
   int *U_prev_row_counts = (int *)calloc(U_n, sizeof(int));

   int *L_done_flags = (int *)calloc(L_nnz, sizeof(int));
   int *U_done_flags = (int *)calloc(U_nnz, sizeof(int));
   int *U_init_flags = (int *)calloc(U_n, sizeof(int));

   int lump = 1;

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = ts->input.num_threads;
      int all_count = 0, L_count = 0, U_init_count = 0, U_count = 0;
      double solve_start;
      int L_nnz_loc = 0;
      int U_nnz_loc = 0;
      int U_n_loc = 0;
      int L_n_loc = 0;
      int kk, ii;

      CSR L_loc;
      CSR U_loc;

      #pragma omp for schedule(static, lump)
      for (int i = 0; i < L_n; i++){
         x[i] = b[i] / L.diag[i];
         y[i] = 0.0;
         L_row_counts[i] = L.i_ptr[i+1] - L.i_ptr[i];
         U_row_counts[i] = U.i_ptr[i+1] - U.i_ptr[i] + 1;
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

      U_loc.i = (int *)calloc(U_nnz_loc, sizeof(int));
      U_loc.j = (int *)calloc(U_nnz_loc, sizeof(int));
      U_loc.data = (double *)calloc(U_nnz_loc, sizeof(double));
      U_loc.diag = (double *)calloc(U_nnz_loc, sizeof(double));
      L_loc.i = (int *)calloc(L_nnz_loc, sizeof(int));
      L_loc.j = (int *)calloc(L_nnz_loc, sizeof(int));
      L_loc.data = (double *)calloc(L_nnz_loc, sizeof(double));
      L_loc.diag = (double *)calloc(L_nnz_loc, sizeof(double));
      double *U_diag_loc = (double *)calloc(U_n_loc, sizeof(double));
      int *L_done_flags_loc = (int *)calloc(L_nnz_loc, sizeof(int));
      int *U_done_flags_loc = (int *)calloc(U_nnz_loc, sizeof(int));
      int *U_init_flags_loc = (int *)calloc(U_n_loc, sizeof(int));

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
  
      solve_start = omp_get_wtime();
      while (1){ 
         if (ts->input.atomic_flag == 1){
            if (ts->input.omp_for_flag == 1){ // TODO: fix bug in this if-statement that causes deadlock
               if (L_count > 0){
                  #pragma omp for nowait schedule(static, lump)
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

                           #pragma omp atomic
                           x[i] -= x[j] * (Lij / Lii);
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
                  #pragma omp for nowait schedule(static, lump)
                  for (int i = 0; i < U_n; i++){
                     int L_row_counts_i;
                     #pragma omp atomic read
                     L_row_counts_i = L_row_counts[i];
                     if (U_init_flags[i] == 0 && L_row_counts_i == 0){
                        double Uii = U.diag[i];

                        #pragma omp atomic
                        y[i] += x[i] / Uii;
                        #pragma omp atomic
                        U_row_counts[i]--;

                        U_init_flags[i] = 1;
                        U_init_count--;
                        all_count--;
                     }
                  }
               }
               if (U_count > 0){
                  #pragma omp for nowait schedule(static, lump)
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

                           #pragma omp atomic
                           y[i] -= y[j] * (Uij / Uii);
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

                           #pragma omp atomic
                           x[i] -= x[j] * (Lij / Lii);
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
                  #pragma omp for nowait schedule(static, lump)
                  for (int i = 0; i < U_n; i++){
                     int L_row_counts_i;
                     #pragma omp atomic read
                     L_row_counts_i = L_row_counts[i];
                     if (U_init_flags_loc[ii] == 0 && L_row_counts_i == 0){
                        double Uii = U_diag_loc[ii];

                        #pragma omp atomic
                        y[i] += x[i] / Uii;
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

                           #pragma omp atomic
                           y[i] -= y[j] * (Uij / Uii);
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
            if (L_count > 0){
               #pragma omp for nowait schedule(static, lump)
               for (int k = 0; k < L_nnz; k++){
                  if (L_done_flags[k] == 0){
                     int i = L.i[k];
                     int j = L.j[k];
                     double Lij = L.data[k];
                     double Lii = L.diag[i];
                     if (L_row_counts[j] == 0){

                        //start = omp_get_wtime();
                        x[i] -= x[j] * (Lij / Lii);
                        #pragma omp atomic
                        L_row_counts[i]--;
                        //ts->output.atomic_wtime[tid] += omp_get_wtime() - start;

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
               #pragma omp for nowait schedule(static, lump)
               for (int i = 0; i < U_n; i++){
                  if (U_init_flags[i] == 0 && L_row_counts[i] == 0){
                     double Uii = U.diag[i];

                     //start = omp_get_wtime();
                     y[i] += x[i] / Uii;
                     #pragma omp atomic
                     U_row_counts[i]--;
                     //ts->output.atomic_wtime[tid] += omp_get_wtime() - start;

                     //dummy_y[i] += dummy_x[i] / Uii;
                     //dummy_U_row_counts[i]--;

                     U_init_flags[i] = 1;
                     U_init_count--;
                     all_count--;
                  }
               }
            }
            if (U_count > 0){
               #pragma omp for nowait schedule(static, lump)
               for (int k = 0; k < U_nnz; k++){
                  if (U_done_flags[k] == 0){
                     int i = U.i[k];
                     int j = U.j[k];
                     double Uij = U.data[k];
                     double Uii = U.diag[i];
                     if (U_row_counts[j] == 0){

                        //start = omp_get_wtime();
                        y[i] -= y[j] * (Uij / Uii);
                        #pragma omp atomic
                        U_row_counts[i]--;
                        //ts->output.atomic_wtime[tid] += omp_get_wtime() - start;

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
   }

   //free((void *)L_row_counts);
   //free((void *)L_prev_row_counts);
   //free((void *)U_row_counts);
   //free((void *)U_prev_row_counts);

   free(L_row_counts);
   free(L_prev_row_counts);
   free(U_row_counts);
   free(U_prev_row_counts);

   free(dummy_L_row_counts);
   free(dummy_U_row_counts);
   free(dummy_x);
   free(dummy_y);

   free(L_done_flags);
   free(U_done_flags);
   free(U_init_flags);

   free(converge_flags);
}
