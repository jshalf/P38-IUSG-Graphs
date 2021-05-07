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

   uint64_t *cycles = (uint64_t *)malloc(num_threads * sizeof(uint64_t));
   double *wtime = (double *)malloc(num_threads * sizeof(double));
   double a = 1.1;
 
   for (int run = 1; run <= num_runs; run++){
      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         double b = (double)tid;

         wtime[tid] = 0;
         cycles[tid] = 0;

         uint64_t cycles_sum = 0;
         #pragma omp barrier

         double wtime_start = omp_get_wtime();
         uint64_t cycles_start = rdtsc();
         if (atomic_flag == 1){
            for (int iters = 0; iters < num_iters; iters++){   
               #pragma omp atomic
               a += b;
            }
         }
         else {
            for (int iters = 0; iters < num_iters; iters++){
               a += b;
            }
         }
         cycles[tid] = rdtsc() - cycles_start;
         wtime[tid] = omp_get_wtime() - wtime_start;
      }

      uint64_t mean_cycles = 0;
      for (int t = 0; t < num_threads; t++){
         mean_cycles += cycles[t];
      }
      mean_cycles /= (uint64_t)num_threads;

      double mean_wtime = 0;
      for (int t = 0; t < num_threads; t++){
         mean_wtime += wtime[t];
      }
      mean_wtime /= (double)num_threads;

      if (verbose_output){
         printf("Wall-clock time %e, num cycles %" PRIu64 "\n", mean_wtime, mean_cycles);
      }
      else {
         printf("%e, %" PRIu64 "\n", mean_wtime, mean_cycles);
      }
   }

   free(cycles);
   free(wtime);

   return 0;
}
