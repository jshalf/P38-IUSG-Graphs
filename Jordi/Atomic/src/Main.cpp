#include "../../src/Main.hpp"
#include "../../src/Misc.hpp"

int main (int argc, char *argv[])
{
   int num_threads = 1;
   int atomic_flag = 1;
   int MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int num_iters = 100;
   int n = 100;

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){
         atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-num_iters") == 0){
         arg_index++;
         num_iters = atoi(argv[arg_index]);
      }
      arg_index++;
   }
   
   omp_set_num_threads(num_threads);

   uint64_t *cycles_sum = (uint64_t *)malloc(num_threads * sizeof(uint64_t));
   double *wtime_sum = (double *)malloc(num_threads * sizeof(double));

   srand(time(NULL));
   double a = RandDouble(-1.0, 1.0);
 
   for (int run = 1; run <= num_runs; run++){
      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         double dummy = 0.0;
         double b = (double)tid;

         uint64_t cycles_sum_loc = 0;
         #pragma omp barrier


         if (atomic_flag == 1){
            uint64_t cycles_start = rdtsc();
            for (int iters = 0; iters < num_iters; iters++){   
               //dummy += b;
               #pragma omp atomic
               a += b;
            }
            uint64_t cycles_stop = rdtsc();
            cycles_sum_loc += cycles_stop - cycles_start;
         }
         else {
            uint64_t cycles_start = rdtsc();
            for (int iters = 0; iters < num_iters; iters++){
               //dummy += b;
               a += b;
            }
            uint64_t cycles_stop = rdtsc();
            cycles_sum_loc += cycles_stop - cycles_start;
         }
         //uint64_t cycles_start = rdtsc();
         //for (int iters = 0; iters < num_iters; iters++){
         //   dummy += b;
         //}
         //uint64_t cycles_stop = rdtsc();
         //cycles_sum_loc -= cycles_stop - cycles_start;


         double wtime_sum_loc = 0;
         #pragma omp barrier


         if (atomic_flag == 1){
            double wtime_start = omp_get_wtime();
            for (int iters = 0; iters < num_iters; iters++){
               //dummy += b;
               #pragma omp atomic
               a += b;
            }
            double wtime_stop = omp_get_wtime();
            wtime_sum_loc += wtime_stop - wtime_start;
         }
         else {
            double wtime_start = omp_get_wtime();
            for (int iters = 0; iters < num_iters; iters++){
               //dummy += b;
               a += b;
            }
            double wtime_stop = omp_get_wtime();
            wtime_sum_loc += wtime_stop - wtime_start;
         }
         //double wtime_start = omp_get_wtime();
         //for (int iters = 0; iters < num_iters; iters++){
         //   dummy += b; 
         //}
         //double wtime_stop = omp_get_wtime();
         //wtime_sum_loc -= wtime_stop - wtime_start;
         
         cycles_sum[tid] = cycles_sum_loc;
         wtime_sum[tid] = wtime_sum_loc;

         //PrintDummy(dummy);
      }

      uint64_t cycles_mean = 0;
      for (int t = 0; t < num_threads; t++){
         cycles_mean += cycles_sum[t];
      }
      cycles_mean /= (uint64_t)num_threads;

      double wtime_mean = 0;
      for (int t = 0; t < num_threads; t++){
         wtime_mean += wtime_sum[t];
      }
      wtime_mean /= (double)num_threads;

      if (verbose_output){
         printf("Wall-clock time  = %e, "
                "num cycles = %" PRIu64 ", "
                "acummulated value = %e\n", 
                wtime_mean,
                cycles_mean,
                a);
      }
      else {
         printf("%e %" PRIu64 " %e\n",
                wtime_mean, cycles_mean, a);
      }
   }

   free(cycles_sum);
   free(wtime_sum);

   return 0;
}
