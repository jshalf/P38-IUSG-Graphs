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

   volatile int *L_row_counts = (int *)calloc(L_n, sizeof(int));
   volatile int *L_prev_row_counts = (int *)calloc(L_n, sizeof(int));
   volatile int *U_row_counts = (int *)calloc(U_n, sizeof(int));
   volatile int *U_prev_row_counts = (int *)calloc(U_n, sizeof(int));

   int *L_done_flags = (int *)calloc(L_nnz, sizeof(int));
   int *U_done_flags = (int *)calloc(U_nnz, sizeof(int));
   int *U_init_flags = (int *)calloc(U_n, sizeof(int));

   int lump = 1;

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = ts->input.num_threads;
      int all_count = 0, L_count = 0, U_init_count = 0, U_count = 0;

      #pragma omp for schedule(static, lump)
      for (int i = 0; i < L_n; i++){
         x[i] = b[i] / L.diag[i];
         y[i] = 0.0;
         L_row_counts[i] = L.i_ptr[i+1] - L.i_ptr[i];
         U_row_counts[i] = U.i_ptr[i+1] - U.i_ptr[i] + 1;
         U_init_count++;
      }
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < L_nnz; k++){
         L_count++;
      }
      #pragma omp for nowait schedule(static, lump)
      for (int k = 0; k < U_nnz; k++){
         U_count++;
      }

      all_count = L_count + U_init_count + U_count;
  
      while (1){ 
         //if (ts->input.atomic_flag == 1){
            if (L_count > 0){
               #pragma omp for nowait schedule(static, lump)
               for (int k = 0; k < L_nnz; k++){
                  if (L_done_flags[k] == 0){
                     int i = L.i[k];
                     int j = L.j[k];
                     double Lij = L.data[k];
                     double Lii = L.diag[i];
                     if (L_row_counts[j] == 0){
                        #pragma atomic
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
                  if (U_init_flags[i] == 0 && L_row_counts[i] == 0){
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
                     int i = U.i[k];
                     int j = U.j[k];
                     double Uij = U.data[k];
                     double Uii = U.diag[i];
                     if (U_row_counts[j] == 0){
                        #pragma atomic
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
         //}
         //else {
         //}

         //if (ts->input.async_flag == 0){
         //   #pragma omp barrier
         //}

         if (all_count == 0){
            break;
         }
      }
   }

   free((void *)L_row_counts);
   free((void *)L_prev_row_counts);
   free((void *)U_row_counts);
   free((void *)U_prev_row_counts);

   free(L_done_flags);
   free(U_done_flags);
   free(U_init_flags);
}
