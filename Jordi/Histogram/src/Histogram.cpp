#include "Histogram.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

int HistogramProgressFunc(Queue *Q,
                          HistogramProgressData *HPD,
                          int flag
                          );

void Histogram_Seq(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   )
{
   for (int i = 0; i < n_index; i++){
      Tally[index[i]]++;
   }
}

void Histogram_Par(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   )
{
   int lump = 1;

   Queue Q, Qp; 
   if (hist->input.MsgQ_flag == 1){
      qAlloc(&Q, n_tally);
      qInitLock(&Q);
      
      int num_threads2 = hist->input.num_threads * hist->input.num_threads;
      qAlloc(&Qp, num_threads2);
      qInitLock(&Qp);
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

      n_index_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         n_index_loc++;
      }
      my_index_part = (int *)calloc(n_index_loc, sizeof(int));
      index_loc = (int *)calloc(n_index_loc, sizeof(int));
      i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         my_index_part[i_loc] = i;
         index_loc[i_loc] = index[i];
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
 
      if (hist->input.MsgQ_flag == 1){
         /*****************
          *   MsgQ wtime
          *****************/
         if (hist->input.MsgQ_wtime_flag == 1){
         }
         /*****************
          *   MsgQ cycles
          *****************/
         else if (hist->input.MsgQ_cycles_flag == 1){
         }
         /*****************
          *  Comp wtime
          *****************/
         else if (hist->input.comp_wtime_flag == 1){
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
            wtime_start = omp_get_wtime();

            for (i_loc = 0; i_loc < n_index_loc; i_loc++){
               int index_i = index_loc[i_loc];
               qPut(&Q, index_i, 1.0);
               num_qPuts++;
               HPD.num_qPuts[index_i % num_threads]++;
            }
            HistogramProgressFunc(&Qp, &HPD, PROGRESS_QPUT);
            int progress_flag;
            do {
               for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
                  int i = my_tally_part[i_loc];
                  double z;
                  int get_flag;
                  get_flag = qGet(&Q, i, &z);
                  if (get_flag == 1){
                     Tally_loc[i_loc]++;
                     num_qGets++;
                     HPD.num_qGets = num_qGets;
                  }
               }
               progress_flag = HistogramProgressFunc(&Qp, &HPD, PROGRESS_QGET);
            } while (progress_flag == 0);

            wtime_stop = omp_get_wtime();
            solve_wtime = wtime_stop - wtime_start;
         }
      }
      else if (hist->input.expand_flag == 1){
         wtime_start = omp_get_wtime();

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }
      else if (hist->input.atomic_flag == 1){
         wtime_start = omp_get_wtime();

         for (i_loc = 0; i_loc < n_index_loc; i_loc++){
            int i = my_index_part[i_loc];
            #pragma omp atomic
            Tally[index[i]]++;
         }

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }
      else {
         wtime_start = omp_get_wtime();

         for (i_loc = 0; i_loc < n_index_loc; i_loc++){
            int i = my_index_part[i_loc];
            Tally[index[i]]++;
         }

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }

      hist->output.solve_wtime_vec[tid] = solve_wtime;

      if (hist->input.MsgQ_flag == 1){
         for (i_loc = 0; i_loc < n_tally_loc; i_loc++){
            int i = my_tally_part[i_loc];
            Tally[i] = Tally_loc[i_loc];
         }

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
      free(index_loc);
   }

   if (hist->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
      
      qDestroyLock(&Qp);
      qFree(&Qp);
   }
}

int HistogramProgressFunc(Queue *Q,
                          HistogramProgressData *HPD,
                          int flag
                          )
{
   int progress_flag = 0;
   int num_threads = HPD->num_threads;
   int tid = HPD->tid;

   if (flag == PROGRESS_QPUT){
      for (int t = 0; t < num_threads; t++){
         //if (t != tid){
            qPut(Q, num_threads*t+tid, HPD->num_qPuts[t]);
         //}
      }
      progress_flag = 1;
   }
   else {
      if (HPD->thread_recv_count < num_threads){
         for (int t = 0; t < num_threads; t++){
            double t_num_gets;
            //if (t != tid){
               int get_flag = qGet(Q, num_threads*tid+t, &t_num_gets);
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
   }

   return progress_flag;
}
