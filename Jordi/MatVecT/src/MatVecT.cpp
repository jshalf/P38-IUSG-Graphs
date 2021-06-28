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
   int num_rows = A.n;
   int lump = 1;

   #pragma omp parallel
   {
      double Axi;

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
   vector<int> row_to_thread;
   int *t_sum, **t_sum_loc;
   /* set up message queues */
   if (mv->input.MsgQ_flag == 1){
      q_size = num_rows;
      qAlloc(&Q, q_size);
      qInitLock(&Q);
      row_to_thread.resize(num_rows);
      t_sum = (int *)calloc(mv->input.num_threads, sizeof(int));
      t_sum_loc = (int **)calloc(mv->input.num_threads, sizeof(int *));
      for (int t = 0; t < mv->input.num_threads; t++){
         t_sum_loc[t] = (int *)calloc(mv->input.num_threads, sizeof(int *));
      }
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = mv->input.num_threads;

      double comp_wtime_start, comp_wtime = 0.0, MsgQ_wtime_start, MsgQ_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles = 0, comp_cycles_start, comp_cycles = 0;
      int num_gets;
      double dummy = 0.0;

      double *y_loc;

      if (mv->input.MsgQ_flag == 1){
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < num_rows; i++){
            row_to_thread[i] = tid;
         }
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < num_rows; i++){
            for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
               int ii = A.j[jj];
               t_sum_loc[tid][row_to_thread[ii]]++;
            }
         }
         #pragma omp for schedule(static, lump)
         for (int t = 0; t < mv->input.num_threads; t++){
            for (int tt = 0; tt < mv->input.num_threads; tt++){
               t_sum[t] += t_sum_loc[tt][t];
            }
         }
         num_gets = t_sum[tid];
         y_loc = (double *)calloc(num_rows, sizeof(double));
      }
  
      double wtime_start = omp_get_wtime();
 
      if (mv->input.MsgQ_flag == 1){
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
                  int j = A.j[jj];
                  MsgQ_wtime_start = omp_get_wtime();
                  qPut(&Q, j, z);
                  MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
               }
            }
            while (num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z, z_accum = 0.0;
                  int get_flag;
                  do {
                     MsgQ_wtime_start = omp_get_wtime();
                     get_flag = qGet(&Q, i, &z);
                     MsgQ_wtime += omp_get_wtime() - MsgQ_wtime_start;
                     if (get_flag == 1){
                        z_accum += z;
                        num_gets--;
                     }
                  } while(get_flag == 1);
                  y_loc[i_loc] += z_accum;
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
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  int j = A.j[jj];
                  MsgQ_cycles_start = rdtsc();
                  qPut(&Q, j, z);
                  MsgQ_cycles += rdtsc() - MsgQ_cycles_start;
               }
            }
            while (num_gets > 0){
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
                        num_gets--;
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
               comp_wtime_start = omp_get_wtime();
               int low = A.start[i];
               int high = A.start[i+1];
               comp_wtime += omp_get_wtime() - comp_wtime_start;
               for (int jj = low; jj < high; jj++){
                  double z;
                  comp_wtime_start = omp_get_wtime();
                  z = A.data[jj] * x[i];
                  int j = A.j[jj];
                  comp_wtime += omp_get_wtime() - comp_wtime_start;
                  qPut(&Q, j, z);
               }
            }
            while (num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z, z_accum = 0.0;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        comp_wtime_start = omp_get_wtime();
                        z_accum += z;
                        num_gets--;
                        comp_wtime += omp_get_wtime() - comp_wtime_start;
                     }
                  } while(get_flag == 1);
                  comp_wtime_start = omp_get_wtime();
                  y_loc[i_loc] += z_accum;
                  i_loc++;
                  comp_wtime += omp_get_wtime() - comp_wtime_start;
               }
            }
         }
         /******************************
          *   MsgQ no-op and comp no-op
          ******************************/
         else if ((mv->input.MsgQ_noop_flag == 1) && (mv->input.comp_noop_flag == 1)){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  int j = A.j[jj];
                  dummy += j;
               }
            }
         }
         /******************************
          *         MsgQ no-op
          ******************************/
         else if (mv->input.MsgQ_noop_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  dummy += z;
               }
            }
            int i_loc = 0;
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               double z, z_accum = (double)rand();
               y_loc[i_loc] += z_accum;
               i_loc++;
            }
         }
         /******************************
          *        Comp no-op
          ******************************/
         else if (mv->input.comp_noop_flag == 1){
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  double z;
                  int j = A.j[jj];
                  qPut(&Q, j, z);
               }
            }
            while (num_gets > 0){
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z, z_accum = 0.0;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        z_accum += z;
                        num_gets--;
                     }
                  } while(get_flag == 1);
               }
            }
         }
         /*************************************************
          * standard scheme (no timers, no no-ops, etc...)
          *************************************************/
         else {
            #pragma omp for schedule(static, lump) nowait
            for (int i = 0; i < num_rows; i++){
               int low = A.start[i];
               int high = A.start[i+1];
               for (int jj = low; jj < high; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  int j = A.j[jj];
                  qPut(&Q, j, z);
               }
            }
            while (num_gets > 0){
               int i_loc = 0;
               #pragma omp for schedule(static, lump) nowait
               for (int i = 0; i < num_rows; i++){
                  double z, z_accum = 0.0;
                  int get_flag;
                  do {
                     get_flag = qGet(&Q, i, &z);
                     if (get_flag == 1){
                        z_accum += z;
                        num_gets--;
                     }
                  } while(get_flag == 1);
                  y_loc[i_loc] += z_accum;
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

         PrintDummy(dummy);
      }
   }

   if (mv->input.MsgQ_flag == 1){
      free(t_sum);
      for (int t = 0; t < mv->input.num_threads; t++){
         free(t_sum_loc[t]);
      }
      free(t_sum_loc);
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
