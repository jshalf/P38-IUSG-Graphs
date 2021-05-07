#include "MatVecT.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

void MatVec_CSR(MatVecData *mv,
                Matrix A, /* sparse matrix */
                double *x, /* vector to be mulitplied with A */
                double *y /* result of Ax */
                )
{
   int lump = 1;

   #pragma omp parallel
   {
      double Axi;
      int num_rows = A.n;

      #pragma omp for schedule(static, lump)
      for (int i = 0; i < num_rows; i++){
         Axi = 0.0;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            Axi += A.data[jj] * x[A.j[jj]];
         }
         y[i] = Axi;
      }
   }
}

void MatVecT_CSR_Seq(MatVecData *mv,
                     Matrix A, /* sparse matrix */
                     double *x, /* vector to be mulitplied with A */
                     double *y /* result of Ax */
                     )
{
   double Axi;
   int num_rows = A.n;

   for (int i = 0; i < num_rows; i++){
      for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
         y[A.j[jj]] += A.data[jj] * x[i];
      }
   }
}

/**************************************************************************
 * Parallel y=A^Tx where y is unknown.
 * A is in compressed sparse row (Matrix) format.
 **************************************************************************/
void MatVecT_CSR(MatVecData *mv,
                 Matrix A, /* sparse matrix */
                 double *x, /* vector to be mulitplied with A */
                 double *y /* result of A^Tx */
                 )
{
   int num_rows = A.n;
   int num_cols = A.m;
   int nnz = A.nnz;

   int lump = 1;

   int q_size;
   Queue Q;
   int *col_counts;
   /* set up message queues */
   if (mv->input.MsgQ_flag == 1){
      q_size = num_rows;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
      col_counts = (int *)calloc(num_rows, sizeof(double));
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = mv->input.num_threads;

      double comp_wtime_start, comp_wtime = 0.0, MsgQ_wtime_start, MsgQ_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles = 0, comp_cycles_start, comp_cycles = 0;

      double *y_loc;
  
      double wtime_start = omp_get_wtime();
 
      if (mv->input.MsgQ_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            int low = A.start[i];
            int high = A.start[i+1];
            for (int jj = low; jj < high; jj++){
               #pragma omp atomic
               col_counts[A.j[jj]]++;
            }
         }
         int my_num_gets = 0;
         int n_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            my_num_gets += col_counts[i];
            n_loc++;
         }
         y_loc = (double *)calloc(num_rows, sizeof(double));

         /*****************
          *   MsgQ wtime
          *****************/
         if (mv->input.MsgQ_wtime_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  MsgQ_wtime_start = omp_get_wtime();
                  qPut(&Q, A.j[jj], z);
                  MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
               }
            }
            while (my_num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  int get_flag;
                  do {
                     MsgQ_wtime_start = omp_get_wtime();
                     get_flag = qGet(&Q, i, &z);
                     MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                     if (get_flag == 1){
                        y_loc[i_loc] += z;
                        my_num_gets--;
                     }
                  } while(get_flag == 1);
                  i_loc++;
               }
            }
         }
         /*****************
          *   MsgQ cycles
          *****************/
         else if (mv->input.MsgQ_cycles_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  MsgQ_cycles_start = rdtsc();
                  qPut(&Q, A.j[jj], z);
                  MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
               }
            }
            while (my_num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  int get_flag;
                  do {
                     MsgQ_cycles_start = rdtsc();
                     get_flag = qGet(&Q, i, &z);
                     MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
                     if (get_flag == 1){
                        y_loc[i_loc] += z;
                        my_num_gets--;
                     }
                  } while(get_flag == 1);
                  i_loc++;
               }
            }
         }
         /*****************
          *  Comp wtime
          *****************/
         else if (mv->input.comp_wtime_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  comp_wtime_start = omp_get_wtime();
                  z = A.data[jj] * x[i];
                  comp_wtime += omp_get_wtime() - comp_wtime_start;
                  qPut(&Q, A.j[jj], z);
               }
            }
            while (my_num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        comp_wtime_start = omp_get_wtime();
                        y_loc[i_loc] += z;
                        comp_wtime += omp_get_wtime() - comp_wtime_start;
                        my_num_gets--;
                     }
                  } while(get_flag == 1);
                  i_loc++;
               }
            }
         }
         /******************************
          *   MsgQ no-op and comp no-op
          ******************************/
         else if ((mv->input.MsgQ_noop_flag == 1) && (mv->input.comp_noop_flag == 1)){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               }
            }
         }
         /******************************
          *         MsgQ no-op
          ******************************/
         else if (mv->input.MsgQ_noop_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
               }
            }
         }
         /******************************
          *        Comp no-op
          ******************************/
         else if (mv->input.comp_noop_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  qPut(&Q, A.j[jj], z);
               }
            }
            while (my_num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        y_loc[i_loc] += z;
                        my_num_gets--;
                     }
                  } while(get_flag == 1);
                  i_loc++;
               }
            }
         }
         /**********************************************
          * standard scheme (no timers, no-ops, etc...)
          **********************************************/
         else {
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  qPut(&Q, A.j[jj], z);
               }
            }
            while (my_num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        y_loc[i_loc] += z;
                        my_num_gets--;
                     }
                  } while(get_flag == 1);
                  i_loc++;
               }
            }
         }
      }
      else if (mv->input.expand_flag == 1){
         int j_offset = num_cols * tid;
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < num_rows; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               mv->y_expand[j_offset + A.j[jj]] += A.data[jj] * x[i];
            }
         }
   
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_cols; i++){
            y[i] = 0;
            for (int j = 0; j < num_threads; j++){
               int jj = j*num_cols + i;
               y[i] += mv->y_expand[jj];
               mv->y_expand[jj] = 0;
            }
         }
      }
      else if (mv->input.atomic_flag == 1){
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               #pragma omp atomic
               y[A.j[jj]] += A.data[jj] * x[i];
            }
         }
      }
      else {
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               y[A.j[jj]] += A.data[jj] * x[i];
            }
         }
      }

      mv->output.solve_wtime_vec[tid] = omp_get_wtime() - wtime_start;

      if (mv->input.MsgQ_flag == 1){
         int i_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            y[i] = y_loc[i_loc];
            i_loc++;
         }
         free(y_loc);

         if (mv->input.MsgQ_wtime_flag == 1){
            mv->output.MsgQ_wtime_vec[tid] = MsgQ_wtime;
         }
         else if (mv->input.MsgQ_cycles_flag == 1){
            mv->output.MsgQ_cycles_vec[tid] = MsgQ_cycles;
         }
         else if (mv->input.comp_wtime_flag == 1){
            mv->output.comp_wtime_vec[tid] = comp_wtime;
         }
         else if (mv->input.comp_cycles_flag == 1){
            mv->output.comp_cycles_vec[tid] = comp_cycles;
         }
      }
   }

   if (mv->input.MsgQ_flag == 1){
      free(col_counts);
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
