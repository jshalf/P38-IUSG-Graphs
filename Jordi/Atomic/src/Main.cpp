#include "../../src/Main.hpp"
#include "../../src/Misc.hpp"

int main (int argc, char *argv[])
{
   int num_threads = 1;
   int atomic_flag = 1;
   int MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int num_iters = 1;
   int n = 100;

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
         n = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-atomic") == 0){
         arg_index++;
         atomic_flag = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){
         MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-num_iters") == 0){
         arg_index++;
         num_iters = atoi(argv[arg_index]);
      }
      arg_index++;
   }
   
   omp_set_num_threads(num_threads);

   //double *x = (double *)calloc(n, sizeof(double));
   //double *b = (double *)calloc(n, sizeof(double));

   uint64_t *cycles = (uint64_t *)malloc(num_threads * sizeof(uint64_t));
   double a = 1.1;
 
   for (int run = 1; run <= num_runs; run++){
      //srand(0);
      //for (int i = 0; i < n; i++){
      //   x[i] = RandDouble(-1.0, 1.0);
      //}
      double start = omp_get_wtime();
      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         double b = (double)tid;
         cycles[tid] = 0;
         if (atomic_flag == 1){
            //uint64_t cycles_start = rdtsc();
            for (int iters = 0; iters < num_iters; iters++){   
               //#pragma omp for nowait
               //for (int i = 0; i < n; i++){
               //   int j = rand() % n;

               //   //double z = RandDouble(-1.0, 1.0);
               //   //#pragma omp atomic
               //   //x[j] += z;

               //   double xj;
               //   #pragma omp atomic read
               //   xj = x[j];
               //   #pragma omp atomic
               //   x[i] += xj;
               //}
               
               uint64_t cycles_start = rdtsc();
               #pragma omp atomic
               a += b;
               uint64_t cycles_stop = rdtsc();
               cycles[tid] += cycles_stop - cycles_start;
            }
            //cycles[tid] = rdtsc() - cycles_start;
         }
         else {
            //uint64_t cycles_start = rdtsc();
            for (int iters = 0; iters < num_iters; iters++){
               //#pragma omp for nowait
               //for (int i = 0; i < n; i++){
               //   int j = rand() % n;

               //   //double z = RandDouble(-1.0, 1.0);
               //   //x[j] += z;

               //   x[i] += x[j];
               //}
               
               uint64_t cycles_start = rdtsc();
               a += b;
               uint64_t cycles_stop = rdtsc();
               cycles[tid] += cycles_stop - cycles_start;
            }
            //cycles[tid] = rdtsc() - cycles_start;
         }
      }
      double wtime = omp_get_wtime() - start;

      uint64_t mean_cycles = 0;
      for (int t = 0; t < num_threads; t++){
         mean_cycles += cycles[t];
      }
      mean_cycles /= (uint64_t)num_threads;

      if (verbose_output){
         printf("Wall-clock time %e, %" PRIu64 "\n", wtime, mean_cycles);
      }
      else {
         printf("%e, %" PRIu64 "\n", wtime, mean_cycles);
      }
   }

   //free(x);
   //free(b);

   free(cycles);

   return 0;
}
