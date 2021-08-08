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

void MatVecT_MsgQ_CSR_Thread(
   Queue *Q,
   Matrix A_loc,
   double *x_loc,
   double *y_loc,
   double *z_loc,
   int *my_rows,
   int n_loc,
   int nnz_loc,
   int num_gets,
   int *num_spins,
   int *num_qPuts,
   int *num_qGets
);

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
      double solve_wtime, wtime_start, wtime_stop;
      double comp_wtime_start, comp_wtime_stop, comp_wtime = 0.0;
      double MsgQ_wtime_start, MsgQ_wtime_stop, MsgQ_wtime = 0.0;
      double MsgQ_put_wtime_start, MsgQ_put_wtime_stop, MsgQ_put_wtime = 0.0;
      double MsgQ_get_wtime_start, MsgQ_get_wtime_stop, MsgQ_get_wtime = 0.0;
      uint64_t MsgQ_cycles_start, MsgQ_cycles_stop, MsgQ_cycles = 0;
      uint64_t MsgQ_put_cycles_start, MsgQ_put_cycles_stop, MsgQ_put_cycles = 0;
      uint64_t MsgQ_get_cycles_start, MsgQ_get_cycles_stop, MsgQ_get_cycles = 0;
      int num_gets, num_spins;
      int num_qGets = 0, num_qPuts = 0;
      int dummy_num_qGets = 0, dummy_num_qPuts = 0;
      int dummy_num_spins = 0;
      int temp_num_qPuts = 0, temp_num_qGets = 0;
      int temp_num_gets = num_gets, temp_num_spins = 0;
      double dummy = 0.0;
      int nnz_loc, n_loc;
      int i_loc, jj_loc;
      double *x_loc, *y_loc, *z_loc;
      int *my_rows;
      Matrix A_loc;

      nnz_loc = 0, n_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < num_rows; i++){
         n_loc++;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            nnz_loc++;
         }
      }
      my_rows = (int *)calloc(n_loc, sizeof(int));
      x_loc = (double *)calloc(n_loc, sizeof(double));
      A_loc.j = (int *)calloc(nnz_loc, sizeof(int));
      A_loc.data = (double *)calloc(nnz_loc, sizeof(double));
      A_loc.start = (int *)calloc(n_loc+1, sizeof(int));
      jj_loc = 0, i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < num_rows; i++){
         my_rows[i_loc] = i;
         x_loc[i_loc] = x[i];
         int i_nnz = 0;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            A_loc.j[jj_loc] = A.j[jj];
            A_loc.data[jj_loc] = A.data[jj];
            jj_loc++;
            i_nnz++;
         }
         A_loc.start[i_loc+1] = A_loc.start[i_loc] + i_nnz;
         i_loc++;
      }

      if (mv->input.MsgQ_flag == 1){
         y_loc = (double *)calloc(n_loc, sizeof(double));
         z_loc = (double *)calloc(nnz_loc, sizeof(double));
         jj_loc = 0;
         i_loc = 0;
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_rows; i++){
            z_loc[i_loc] = RandDouble(0.0, 1.0);
            i_loc++;
         }

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
      }
      else if (mv->input.expand_flag == 1){
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < num_cols; i++){
            for (int j = 0; j < num_threads; j++){
               int jj = j*num_cols + i;
               mv->y_expand[jj] = 0;
            }
         }
      }

      #pragma omp barrier
 
      if (mv->input.MsgQ_flag == 1){
         /*****************
          *   MsgQ wtime
          *****************/
         if (mv->input.MsgQ_wtime_flag == 1){
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }
            MsgQ_wtime_start = omp_get_wtime();
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               qPut(&Q, j, z);
               num_qPuts++;
            }
            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_put_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;
            num_spins = 0;
            temp_num_gets = num_gets;
            MsgQ_wtime_start = omp_get_wtime();
            while (temp_num_gets > 0){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z, z_accum = 0.0;
                  int get_flag;
                  get_flag = qGet(&Q, i, &z);
                  if (get_flag == 1){
                     num_qGets++;
                     z_accum += z;
                     temp_num_gets--;
                  }
                  y_loc[i_loc] += z_accum;
               }
               num_spins++;
            }
            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_get_wtime += MsgQ_wtime_stop - MsgQ_wtime_start;

            /* re-do some computation to re-create cache behaviour */
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }

            /* now compute wtime when MsgQs are removed */
            MsgQ_wtime_start = omp_get_wtime();
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               dummy += (z + j);
               dummy_num_qPuts++;
            }
            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_put_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;
            dummy_num_spins = 0;
            MsgQ_wtime_start = omp_get_wtime();
            for (int s = 0; s < num_spins; s++){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z = 1.1, z_accum = 1.1;
                  int get_flag;
                  dummy_num_qGets++;
                  z_accum += z;
                  y_loc[i_loc] += z_accum;
               }
               dummy_num_spins++;
            }
            /* subtract elapsed MsgQ-less wtime */
            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_get_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;
         }
         /*****************
          *   MsgQ cycles
          *****************/
         else if (mv->input.MsgQ_cycles_flag == 1){
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }
            MsgQ_cycles_start = rdtsc();
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               qPut(&Q, j, z);
               num_qPuts++;
            }
            MsgQ_cycles_stop = rdtsc();
            MsgQ_put_cycles += MsgQ_cycles_stop - MsgQ_cycles_start;
            num_spins = 0;
            temp_num_gets = num_gets;
            MsgQ_cycles_start = rdtsc();
            while (temp_num_gets > 0){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z, z_accum = 0.0;
                  int get_flag;
                  get_flag = qGet(&Q, i, &z);
                  if (get_flag == 1){
                     num_qGets++;
                     z_accum += z;
                     temp_num_gets--;
                  }
                  y_loc[i_loc] += z_accum;
               }
               num_spins++;
            }
            MsgQ_cycles_stop = rdtsc();
            MsgQ_get_cycles += MsgQ_cycles_stop - MsgQ_cycles_start;

            /* re-do some computation to re-create cache behaviour */
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }

            /* now compute cycles when MsgQs are removed */
            MsgQ_cycles_start = rdtsc();
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               dummy += (z + j);
               dummy_num_qPuts++;
            }
            MsgQ_cycles_stop = rdtsc();
            MsgQ_put_cycles -= MsgQ_cycles_stop - MsgQ_cycles_start;
            dummy_num_spins = 0;
            MsgQ_cycles_start = rdtsc();
            for (int s = 0; s < num_spins; s++){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z = 1.1, z_accum = 1.1;
                  int get_flag;
                  dummy_num_qGets++;
                  z_accum += z;
                  y_loc[i_loc] += z_accum;
               }
               dummy_num_spins++;
            }
            /* subtract elapsed MsgQ-less cycles */
            MsgQ_cycles_stop = rdtsc();
            MsgQ_get_cycles -= MsgQ_cycles_stop - MsgQ_cycles_start;
         }
         /*****************
          *  Comp wtime
          *****************/
         else if (mv->input.comp_wtime_flag == 1){
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               qPut(&Q, j, z);
               num_qPuts++;
            }
            num_spins = 0;
            temp_num_gets = num_gets;
            while (temp_num_gets > 0){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z, z_accum = 0.0;
                  int get_flag;
                  get_flag = qGet(&Q, i, &z);
                  if (get_flag == 1){
                     num_qGets++;
                     z_accum += z;
                     temp_num_gets--;
                  }
                  y_loc[i_loc] += z_accum;
               }
               num_spins++;
            }

            comp_wtime_start = omp_get_wtime();
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               for (jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
                  z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
               }
            }
            for (jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
               int j = A_loc.j[jj_loc];
               double z = z_loc[jj_loc];
               dummy += (z + j);
               dummy_num_qPuts++;
            }
            dummy_num_spins = 0;
            for (int s = 0; s < num_spins; s++){
               for (int i_loc = 0; i_loc < n_loc; i_loc++){
                  int i = my_rows[i_loc];
                  double z = 1.1, z_accum = 1.1;
                  dummy_num_qGets++;
                  z_accum += z;
                  y_loc[i_loc] += z_accum;
               }
               dummy_num_spins++;
            }
            comp_wtime_stop = omp_get_wtime();
            comp_wtime += comp_wtime_stop - comp_wtime_start;
         }
         /******************************
          *   MsgQ no-op and comp no-op
          ******************************/
         else if ((mv->input.MsgQ_noop_flag == 1) && (mv->input.comp_noop_flag == 1)){
         }
         /******************************
          *         MsgQ no-op
          ******************************/
         else if (mv->input.MsgQ_noop_flag == 1){
         }
         /******************************
          *        Comp no-op
          ******************************/
         else if (mv->input.comp_noop_flag == 1){
         }
         /*************************************************
          * standard scheme (no timers, no no-ops, etc...)
          *************************************************/
         else {
            //MatVecT_MsgQ_CSR_Thread(
            //   &Q, A_loc, x_loc, y_loc, z_loc, my_rows, n_loc, nnz_loc,
            //   temp_num_gets, &temp_num_spins, &temp_num_qPuts, &temp_num_qGets);
            //#pragma omp barrier

            wtime_start = omp_get_wtime();
            MatVecT_MsgQ_CSR_Thread(
               &Q, A_loc, x_loc, y_loc, z_loc, my_rows, n_loc, nnz_loc,
               num_gets, &num_spins, &num_qPuts, &num_qGets);
            wtime_stop = omp_get_wtime();
            solve_wtime = wtime_stop - wtime_start;
         }
      }
      else if (mv->input.expand_flag == 1){
         wtime_start = omp_get_wtime();
         int j_offset = num_cols * tid;

         //#pragma omp for schedule(static, lump)
         //for (int i = 0; i < num_rows; i++){
         //   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
         //      mv->y_expand[j_offset + A.j[jj]] += A.data[jj] * x[i];
         //   }
         //}

         jj_loc = 0;
         for (int i_loc = 0; i_loc < n_loc; i_loc++){
            for (int jj = A_loc.start[i_loc]; jj < A_loc.start[i_loc+1]; jj++){
               mv->y_expand[j_offset + A_loc.j[jj_loc]] += A_loc.data[jj_loc] * x_loc[i_loc];
               jj_loc++;
            }
         }
         #pragma omp barrier
   
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < num_cols; i++){
            y[i] = 0;
            for (int j = 0; j < num_threads; j++){
               int jj = j*num_cols + i;
               y[i] += mv->y_expand[jj];
            }
         }
         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }
      else if (mv->input.atomic_flag == 1){
         wtime_start = omp_get_wtime();

         //#pragma omp for schedule(static, lump) nowait
         //for (int i = 0; i < num_rows; i++){
         //   for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
         //      #pragma omp atomic
         //      y[A.j[jj]] += A.data[jj] * x[i];
         //   }
         //}

         jj_loc = 0;
         for (int i_loc = 0; i_loc < n_loc; i_loc++){
            for (int jj = A_loc.start[i_loc]; jj < A_loc.start[i_loc+1]; jj++){
               #pragma omp atomic
               y[A_loc.j[jj_loc]] += A_loc.data[jj_loc] * x_loc[i_loc];
               jj_loc++;
            }
         }
         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }
      else {
         wtime_start = omp_get_wtime();

        // #pragma omp for schedule(static, lump) nowait
        // for (int i = 0; i < num_rows; i++){
        //    for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
        //       y[A.j[jj]] += A.data[jj] * x[i];
        //    }
        // }

         jj_loc = 0;
         for (int i_loc = 0; i_loc < n_loc; i_loc++){
            for (int jj = A_loc.start[i_loc]; jj < A_loc.start[i_loc+1]; jj++){
               y[A_loc.j[jj_loc]] += A_loc.data[jj_loc] * x_loc[i_loc];
               jj_loc++;
            }
         }

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }

      mv->output.solve_wtime_vec[tid] = solve_wtime;

      if (mv->input.MsgQ_flag == 1){
         for (int i_loc = 0; i_loc < n_loc; i_loc++){
            int i = my_rows[i_loc];
            y[i] = y_loc[i_loc];
         }
         free(y_loc);
         free(z_loc);

         //if (mv->input.MsgQ_wtime_flag == 1){
            mv->output.MsgQ_put_wtime_vec[tid] = MsgQ_put_wtime;
            mv->output.MsgQ_get_wtime_vec[tid] = MsgQ_get_wtime;
         //}
         //else if (mv->input.MsgQ_cycles_flag == 1){
            mv->output.MsgQ_put_cycles_vec[tid] = MsgQ_put_cycles;
            mv->output.MsgQ_get_cycles_vec[tid] = MsgQ_get_cycles;
         //}
         //else if (mv->input.comp_wtime_flag == 1){
            mv->output.comp_wtime_vec[tid] = comp_wtime;
         //}

         mv->output.num_qGets_vec[tid] = num_qGets;
         mv->output.num_qPuts_vec[tid] = num_qPuts;

         dummy += (dummy_num_spins + dummy_num_qPuts + dummy_num_qGets);
         PrintDummy(dummy);

         #pragma omp barrier
      }

      free(my_rows);
      free(x_loc);
      free(A_loc.j);
      free(A_loc.data);
      free(A_loc.start);
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

void MatVecT_MsgQ_CSR_Thread(
   Queue *Q,
   Matrix A_loc,
   double *x_loc,
   double *y_loc,
   double *z_loc,
   int *my_rows,
   int n_loc,
   int nnz_loc,
   int num_gets,
   int *num_spins,
   int *num_qPuts,
   int *num_qGets)
{
   for (int i_loc = 0; i_loc < n_loc; i_loc++){
      for (int jj_loc = A_loc.start[i_loc]; jj_loc < A_loc.start[i_loc+1]; jj_loc++){
         z_loc[jj_loc] = A_loc.data[jj_loc] * x_loc[i_loc];
      }
   }
   for (int jj_loc = 0; jj_loc < nnz_loc; jj_loc++){
      int j = A_loc.j[jj_loc];
      double z = z_loc[jj_loc];
      qPut(Q, j, z);
      num_qPuts[0]++;
   }
   num_spins[0] = 0;
   while (num_gets > 0){
      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         int i = my_rows[i_loc];
         double z, z_accum = 0.0;
         int get_flag;
         get_flag = qGet(Q, i, &z);
         if (get_flag == 1){
            num_qGets[0]++;
            z_accum += z;
            num_gets--;
         }
         y_loc[i_loc] += z_accum;
      }
      num_spins[0]++;
   }
}
