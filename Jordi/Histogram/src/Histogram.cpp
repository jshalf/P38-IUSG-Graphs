#include "Histogram.hpp"
#include "Matrix.hpp"
#include "MsgQ.hpp"
#include "Misc.hpp"

int HistogramProgressFunc(MessageQueue<int> *Q,
                          HistogramProgressData *HPD,
                          int flag
                          );
int HistogramProgressFunc_Comp(MessageQueue<int> *Q,
                               HistogramProgressData *HPD,
                               int flag,
                               double *dummy
                               );


void Histogram_Seq(HistogramData *hist,
                   vector<int> index, int n_index,
                   int *Tally, int n_tally
                   )
{
   for (int i = 0; i < n_index; i++){
      Tally[index[i]]++;
   }
}

void Histogram_Par(HistogramData *hist,
                   vector<int> index, int n_index,
                   int *Tally, int n_tally
                   )
{
   int lump = 1;

   MessageQueue<int> *Q, *Qp; 
   if (hist->input.MsgQ_flag == 1){
      int num_threads2 = hist->input.num_threads * hist->input.num_threads;

      Q = new MessageQueue<int>(n_tally);
      Qp = new MessageQueue<int>(num_threads2);
   }

   double *Tally_reduce;
   if (hist->input.reduce_flag == 1){
      Tally_reduce = (double *)calloc(hist->input.num_threads * n_tally, sizeof(double));
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = hist->input.num_threads;
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
      int n_index_loc, n_tally_loc;
      int i_loc, jj_loc;
      int *my_index_part, *my_tally_part;
      int *index_loc, *Tally_loc;
      HistogramProgressData HPD;

      wtime_start = omp_get_wtime();
      n_index_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         n_index_loc++;
      }
      my_index_part = (int *)calloc(n_index_loc, sizeof(int));
      //index_loc = (int *)calloc(n_index_loc, sizeof(int));
      i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         my_index_part[i_loc] = i;
         //index_loc[i_loc] = index[i];
         i_loc++;
      }

      n_tally_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_tally; i++){
         n_tally_loc++;
      }
      my_tally_part = (int *)calloc(n_tally_loc, sizeof(int));
      Tally_loc = (int *)calloc(n_tally_loc, sizeof(int));
      i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_tally; i++){
         my_tally_part[i_loc] = i;
         i_loc++;
      }

      if (hist->input.MsgQ_flag == 1){
         HPD.thread_recv_count = 0;
         HPD.actual_num_qGets = 0;
         HPD.num_qGets = 0;
         HPD.num_qPuts = (int *)calloc(num_threads, sizeof(int));
         HPD.tid = tid;
         HPD.num_threads = num_threads;
      }

      #pragma omp barrier

      wtime_start = omp_get_wtime(); 
      if (hist->input.MsgQ_flag == 1){
         /*****************
          *   MsgQ wtime
          *****************/
         if (hist->input.MsgQ_wtime_flag == 1){
            MsgQ_wtime_start = omp_get_wtime();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               Q->qPut(index_i, 1);
               num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc(Qp, &HPD, PROGRESS_QPUT);
            int num_spins = 0;
            int progress_flag;
            do {
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  int z;
                  int get_flag;
                  get_flag = Q->qGet(i, &z);
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc(Qp, &HPD, PROGRESS_QGET);
               num_spins++;
            } while (progress_flag == 0);
            num_qGets += num_threads;
            num_qPuts += num_threads;

            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_wtime = MsgQ_wtime_stop - MsgQ_wtime_start;


            MsgQ_wtime_start = omp_get_wtime();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               dummy += index_i;
               dummy_num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QPUT, &dummy);
            for (int s = 0; s < num_spins; s++){
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  dummy += i;
                  double z;
                  int get_flag;
                  get_flag = 1;
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QGET, &dummy);
               dummy_num_spins++;
            }
            num_qGets += num_threads;
            num_qPuts += num_threads;

            MsgQ_wtime_stop = omp_get_wtime();
            MsgQ_wtime -= MsgQ_wtime_stop - MsgQ_wtime_start;
         }
         /*****************
          *   MsgQ cycles
          *****************/
         else if (hist->input.MsgQ_cycles_flag == 1){
            MsgQ_cycles_start = rdtsc();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               Q->qPut(index_i, 1);
               num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc(Qp, &HPD, PROGRESS_QPUT);
            int num_spins = 0;
            int progress_flag;
            do {
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  int z;
                  int get_flag;
                  get_flag = Q->qGet(i, &z);
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc(Qp, &HPD, PROGRESS_QGET);
               num_spins++;
            } while (progress_flag == 0);
            num_qGets += num_threads;
            num_qPuts += num_threads;

            MsgQ_cycles_stop = rdtsc();
            MsgQ_cycles = MsgQ_cycles_stop - MsgQ_cycles_start;


            MsgQ_cycles_start = rdtsc();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               dummy += index_i;
               dummy_num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QPUT, &dummy);
            for (int s = 0; s < num_spins; s++){
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  dummy += i;
                  int z;
                  int get_flag;
                  get_flag = 1;
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QGET, &dummy);
               dummy_num_spins++;
            }
            num_qGets += num_threads;
            num_qPuts += num_threads;

            MsgQ_cycles_stop = rdtsc();
            MsgQ_cycles -= MsgQ_cycles_stop - MsgQ_cycles_start;
         }
         /*****************
          *  Comp wtime
          *****************/
         else if (hist->input.comp_wtime_flag == 1){
            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               Q->qPut(index_i, 1);
               num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc(Qp, &HPD, PROGRESS_QPUT);
            int num_spins = 0;
            int progress_flag;
            do {
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  int z;
                  int get_flag;
                  get_flag = Q->qGet(i, &z);
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc(Qp, &HPD, PROGRESS_QGET);
               num_spins++;
            } while (progress_flag == 0);
            num_qGets += num_threads;
            num_qPuts += num_threads;


            comp_wtime_start = omp_get_wtime();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc];
               int index_i = index[i];
               dummy += index_i;
               dummy_num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QPUT, &dummy);
            for (int s = 0; s < num_spins; s++){
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  dummy += i;
                  double z;
                  int get_flag;
                  get_flag = 1;
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     dummy_num_qGets++;
                     HPD.num_qGets = dummy_num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc_Comp(Qp, &HPD, PROGRESS_QGET, &dummy);
               dummy_num_spins++;
            }
            num_qGets += num_threads;
            num_qPuts += num_threads;

            comp_wtime_stop = omp_get_wtime();
            comp_wtime = comp_wtime_stop - comp_wtime_start;
         }
         /******************************
          *   MsgQ no-op and comp no-op
          ******************************/
         else if ((hist->input.MsgQ_noop_flag == 1) && (hist->input.comp_noop_flag == 1)){
         }
         /******************************
          *         MsgQ no-op
          ******************************/
         else if (hist->input.MsgQ_noop_flag == 1){
         }
         /******************************
          *        Comp no-op
          ******************************/
         else if (hist->input.comp_noop_flag == 1){
         }
         /*************************************************
          * standard scheme (no timers, no no-ops, etc...)
          *************************************************/
         else {
            /* first send indices to be tallied by other threads */
            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int i = my_index_part[i_loc]; 
               int index_i = index[i];
               Q->qPut(index_i, 1);
               num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            /* initiate non-blocking progress */
            HistogramProgressFunc(Qp, &HPD, PROGRESS_QPUT);
            int progress_flag;
            int num_spins = 0;
            do {
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  int z;
                  int get_flag;
                  /* if data arrives, update local tally */
                  get_flag = Q->qGet(i, &z);
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               /* check and advance non-blocking progress */
               progress_flag = HistogramProgressFunc(Qp, &HPD, PROGRESS_QGET);
            } while (progress_flag == 0);
            num_qGets += num_threads;
            num_qPuts += num_threads;
         }
      }
      else if (hist->input.reduce_flag == 1){
         int j_offset = n_tally * tid;
         /* update local tally */
         #pragma omp for schedule(static, lump)
         for (int i = 0; i < n_index; i++){
            Tally_reduce[j_offset + index[i]]++;
         }
         /* sum local tallys */
         #pragma omp for schedule(static, lump) nowait
         for (int i = 0; i < n_tally; i++){
            Tally[i] = 0;
            for (int j = 0; j < num_threads; j++){
               int jj = j*n_tally + i;
               Tally[i] += Tally_reduce[jj];
            }
         }
      }
      else if (hist->input.atomic_flag == 1){
         /* atomically update global tally */
         for (i_loc = 0; i_loc < n_index_loc; i_loc++){
            int i = my_index_part[i_loc];
            #pragma omp atomic
            Tally[index[i]]++;
         }
      }
      else {
         for (i_loc = 0; i_loc < n_index_loc; i_loc++){
            int i = my_index_part[i_loc];
            Tally[index[i]]++;
         }
      }

      wtime_stop = omp_get_wtime();
      solve_wtime = wtime_stop - wtime_start;

      hist->output.solve_wtime_vec[tid] = solve_wtime;

      if (hist->input.MsgQ_flag == 1){
         for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
            int i = my_tally_part[i_loc];
            Tally[i] = Tally_loc[i_loc];
         }

         MsgQ_put_wtime = MsgQ_wtime;

         //if (hist->input.MsgQ_wtime_flag == 1){
            hist->output.MsgQ_put_wtime_vec[tid] = MsgQ_put_wtime;
            hist->output.MsgQ_get_wtime_vec[tid] = MsgQ_get_wtime;
         //}
         //else if (hist->input.MsgQ_cycles_flag == 1){
            hist->output.MsgQ_put_cycles_vec[tid] = MsgQ_put_cycles;
            hist->output.MsgQ_get_cycles_vec[tid] = MsgQ_get_cycles;
         //}
         //else if (hist->input.comp_wtime_flag == 1){
            hist->output.comp_wtime_vec[tid] = comp_wtime;
         //}

         hist->output.num_qGets_vec[tid] = num_qGets;
         hist->output.num_qPuts_vec[tid] = num_qPuts;

         dummy += (dummy_num_spins + dummy_num_qPuts + dummy_num_qGets);
         PrintDummy(dummy);

         free(HPD.num_qPuts);
         free(Tally_loc);
         #pragma omp barrier
      }

      free(my_index_part);
      free(my_tally_part);
      //free(index_loc);
   }

   if (hist->input.reduce_flag == 1){
      free(Tally_reduce);
   }

   if (hist->input.MsgQ_flag == 1){
   }
}

/* progress func */
int HistogramProgressFunc(MessageQueue<int> *Q,
                          HistogramProgressData *HPD,
                          int flag
                          )
{
   int progress_flag = 0;
   int num_threads = HPD->num_threads;
   int tid = HPD->tid;

   /* initial call should have the PROGRESS_QPUT flag */
   if (flag == PROGRESS_QPUT){
      for (int t = 0; t < num_threads; t++){
         //if (t != tid){
            Q->qPut(num_threads*t+tid, HPD->num_qPuts[t]);
         //}
      }
      progress_flag = 1;
   }
   /* subsequent calls should have the PROGRESS_QGET flag */
   else {
      if (HPD->thread_recv_count < num_threads){
         /* check for progress messages from other threads */
         for (int t = 0; t < num_threads; t++){
            int t_num_gets;
            //if (t != tid){
               int get_flag = Q->qGet(num_threads*tid+t, &t_num_gets);
               if (get_flag == 1){
                  /* upon receipt, update total number of messages that
                     need to be received to complete the histogram computation */
                  HPD->actual_num_qGets += (int)t_num_gets;
                  HPD->thread_recv_count++;
               }
            //}
         }
      }
      /* if we've received from all threads, no need to keep checking the progress queue */
      if (HPD->thread_recv_count == num_threads){
         if (HPD->num_qGets == HPD->actual_num_qGets){
            progress_flag = 1;
         }
      }
   }

   return progress_flag;
}

/* progress func where only the computation is performed */
int HistogramProgressFunc_Comp(MessageQueue<int> *Q,
                               HistogramProgressData *HPD,
                               int flag,
                               double *dummy 
                               )
{
   int progress_flag = 0;
   int num_threads = HPD->num_threads;
   int tid = HPD->tid;
   double dummy_temp = 0.0;

   if (flag == PROGRESS_QPUT){
      for (int t = 0; t < num_threads; t++){
         //if (t != tid){
            *dummy += (num_threads*t+tid + HPD->num_qPuts[t]);
         //}
      }
      progress_flag = 1;
   }
   else {
      if (HPD->thread_recv_count < num_threads){
         for (int t = 0; t < num_threads; t++){
            double t_num_gets;
            //if (t != tid){
               int get_flag = 1;//qGet(Q, num_threads*tid+t, &t_num_gets);
               if (get_flag == 1){
                  HPD->actual_num_qGets += (int)t_num_gets;
                  HPD->thread_recv_count++;
               }
            //}
         }
      }
      if (HPD->thread_recv_count == num_threads){
         if (HPD->num_qGets == HPD->actual_num_qGets){
            progress_flag = 1;
         }
      }
      *dummy += progress_flag + HPD->thread_recv_count + HPD->actual_num_qGets + HPD->num_qGets;
   }

   return progress_flag;
}
