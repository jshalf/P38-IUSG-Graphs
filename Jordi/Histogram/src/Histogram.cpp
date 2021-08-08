#include "Histogram.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"
#include "../../src/Misc.hpp"

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
      int *my_index;

      n_index_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         n_index_loc++;
      }
      n_tally_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_tally; i++){
         n_tally_loc++;
      }
      my_index = (int *)calloc(n_index_loc, sizeof(int));
      i_loc = 0;
      #pragma omp for schedule(static, lump) nowait
      for (int i = 0; i < n_index; i++){
         my_index[i_loc] = i;
         i_loc++;
      }

      if (hist->input.MsgQ_flag == 1){
         // compute num_gets here
      }
      else if (hist->input.expand_flag == 1){
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
            int i = my_index[i_loc];
            #pragma omp atomic
            Tally[index[i]]++;
         }

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }
      else {
         wtime_start = omp_get_wtime();

         for (i_loc = 0; i_loc < n_index_loc; i_loc++){
            int i = my_index[i_loc];
            Tally[index[i]]++;
         }

         wtime_stop = omp_get_wtime();
         solve_wtime = wtime_stop - wtime_start;
      }

      hist->output.solve_wtime_vec[tid] = solve_wtime;

      if (hist->input.MsgQ_flag == 1){
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

         #pragma omp barrier
      }

      free(my_index);
   }

   if (hist->input.MsgQ_flag == 1){
   }
}
