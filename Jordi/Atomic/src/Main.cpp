#include "../../src/Main.hpp"
#include "../../src/Misc.hpp"

#define ATOMIC_MULTIPLE_ACCUM 0
#define ATOMIC_SINGLE_WRITE_MULTIPLE_READ 1

int main (int argc, char *argv[])
{
   int num_threads = 1;
   int atomic_flag = 1;
   int MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int num_iters = 100;
   int n = 0;
   int atomic_scheme = ATOMIC_MULTIPLE_ACCUM;

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
      else if (strcmp(argv[arg_index], "-atomic_scheme") == 0){
         arg_index++;
         if (strcmp(argv[arg_index], "ma") == 0){ /* multiple accumulate */
            atomic_scheme = ATOMIC_MULTIPLE_ACCUM;
         }
         else if (strcmp(argv[arg_index], "swmr") == 0){ /* single write multiple read */
            atomic_scheme = ATOMIC_SINGLE_WRITE_MULTIPLE_READ;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
         n = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-atomic_scheme <string>:        atomic scheme.\n");
      printf("    multiple accumulate:        level scheduled method.\n");
      printf("    single-write multiple-read: asynchronous method.\n");
      printf("-n <int value>:                 size of array.\n");
      printf("-num_threads <int value>:       number of OpenMP threads.\n");
      printf("-no_atomic:                     turn off atomics.  Only meant for performance measurements and will not produce a correct result.\n");
      printf("-num_runs <int value>:          number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                      verbose output.\n");
      printf("-num_iters <int value>:         number of iterations in which each thread performs an atomic\n");
      return 0;
   }

   int num_threads_write, num_threads_read;
   if (atomic_scheme == ATOMIC_SINGLE_WRITE_MULTIPLE_READ){
      num_threads_write = 1;
      num_threads_read = num_threads-num_threads_write;
   }
   else {
      num_threads_write = num_threads;
      num_threads_read = num_threads;
   }

   srand(time(NULL));
   int root_tid = 0;//RandInt(0, num_threads-1, 0.0);
   //printf("root tid %d\n", root_tid);
   
   omp_set_num_threads(num_threads);

   uint64_t *read_cycles_sum = (uint64_t *)calloc(num_threads, sizeof(uint64_t));
   double *read_wtime_sum = (double *)calloc(num_threads, sizeof(double));
   uint64_t *write_cycles_sum = (uint64_t *)calloc(num_threads, sizeof(uint64_t));
   double *write_wtime_sum = (double *)calloc(num_threads, sizeof(double));
   uint64_t *accum_cycles_sum = (uint64_t *)calloc(num_threads, sizeof(uint64_t));
   double *accum_wtime_sum = (double *)calloc(num_threads, sizeof(double));

   srand(time(NULL));
   int m = max(1, n * num_threads);
   double *a = (double *)calloc(m, sizeof(double));
   // first touch
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      for (int i = 0; i < n; i++){ 
         a[n*tid+i] = RandDouble(-1.0, 1.0);
      }
   }
 
   for (int run = 1; run <= num_runs; run++){
      //int root_tid = RandInt(0, num_threads-1, 0.0);
      #pragma omp parallel
      {
         int tid = omp_get_thread_num();
         double dummy = 0.0;
         double b = (double)tid;

         uint64_t cycles_sum_loc = 0;
         #pragma omp barrier


         if (atomic_scheme == ATOMIC_SINGLE_WRITE_MULTIPLE_READ){
            /* First measure the number of cycles */
            if (atomic_flag == 1){
               /* Root thread writes to array */
               if (tid == root_tid){
                  uint64_t cycles_start = rdtsc();
                  for (int iters = 0; iters < num_iters; iters++){
                     #pragma omp atomic write relaxed
                     a[n*tid] = iters;
                  }
                  uint64_t cycles_stop = rdtsc();
                  cycles_sum_loc += cycles_stop - cycles_start;
               }
               /* All other threads read from array */
               else {
                  uint64_t cycles_start = rdtsc();
                  for (int iters = 0; iters < num_iters; iters++){
                     #pragma omp atomic read relaxed
                     b = a[n*tid];
                  }
                  uint64_t cycles_stop = rdtsc();
                  cycles_sum_loc += cycles_stop - cycles_start;
               }
            }
            else {
               /* Root thread writes to array */
               if (tid == root_tid){
                  uint64_t cycles_start = rdtsc();
                  for (int iters = 0; iters < num_iters; iters++){
                     //#pragma omp atomic write
                     a[n*tid] = iters;
                  }
                  uint64_t cycles_stop = rdtsc();
                  cycles_sum_loc += cycles_stop - cycles_start;
               }
               /* All other threads read from array */
               else {
                  uint64_t cycles_start = rdtsc();
                  for (int iters = 0; iters < num_iters; iters++){
                     //#pragma omp atomic read
                     b = a[n*tid];
                  }
                  uint64_t cycles_stop = rdtsc();
                  cycles_sum_loc += cycles_stop - cycles_start;
               }
            }

            double wtime_sum_loc = 0;
            #pragma omp barrier

            /* Now measure wall-clock time */
            if (atomic_flag == 1){
               /* Root thread writes to array */
               if (tid == root_tid){
                  double wtime_start = omp_get_wtime();
                  for (int iters = 0; iters < num_iters; iters++){
                     #pragma omp atomic write relaxed
                     a[n*tid] = b;
                  }
                  double wtime_stop = omp_get_wtime();
                  wtime_sum_loc += wtime_stop - wtime_start;
               }
               /* All other threads read from array */
               else {
                  double wtime_start = omp_get_wtime();
                  for (int iters = 0; iters < num_iters; iters++){
                     #pragma omp atomic read relaxed
                     b = a[n*tid];
                  }
                  double wtime_stop = omp_get_wtime();
                  wtime_sum_loc += wtime_stop - wtime_start;
               }
            }
            else {
               /* Root thread writes to array */
               if (tid == root_tid){
                  double wtime_start = omp_get_wtime();
                  for (int iters = 0; iters < num_iters; iters++){
                     //#pragma omp atomic write
                     a[n*tid] = b;
                  }
                  double wtime_stop = omp_get_wtime();
                  wtime_sum_loc += wtime_stop - wtime_start;
               }
               /* All other threads read from array */
               else {
                  double wtime_start = omp_get_wtime();
                  for (int iters = 0; iters < num_iters; iters++){
                     //#pragma omp atomic read
                     b = a[n*tid];
                  }
                  double wtime_stop = omp_get_wtime();
                  wtime_sum_loc += wtime_stop - wtime_start;
               }
            }

            if (tid == root_tid){
               write_cycles_sum[tid] = cycles_sum_loc;
               write_wtime_sum[tid] = wtime_sum_loc;
            }
            else {
               read_cycles_sum[tid] = cycles_sum_loc;
               read_wtime_sum[tid] = wtime_sum_loc;
            }

            #pragma omp atomic
            a[n*tid] += b;
         }
         else {
            /* First measure the number of cycles */
            if (atomic_flag == 1){
               uint64_t cycles_start = rdtsc();
               for (int iters = 0; iters < num_iters; iters++){   
                  //dummy += b;
                  /* All threads atomically update array */
                  #pragma omp atomic
                  a[n*tid] += b;
               }
               uint64_t cycles_stop = rdtsc();
               cycles_sum_loc += cycles_stop - cycles_start;
            }
            else {
               uint64_t cycles_start = rdtsc();
               for (int iters = 0; iters < num_iters; iters++){
                  //dummy += b;
                  /* All threads atomically update array */
                  a[n*tid] += b;
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

            /* First measure the number of cycles */
            if (atomic_flag == 1){
               double wtime_start = omp_get_wtime();
               for (int iters = 0; iters < num_iters; iters++){
                  //dummy += b;
                  /* All threads atomically update array */
                  #pragma omp atomic
                  a[n*tid] += b;
               }
               double wtime_stop = omp_get_wtime();
               wtime_sum_loc += wtime_stop - wtime_start;
            }
            else {
               double wtime_start = omp_get_wtime();
               for (int iters = 0; iters < num_iters; iters++){
                  //dummy += b;
                  /* All threads atomically update array */
                  a[n*tid] += b;
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
            
            accum_cycles_sum[tid] = cycles_sum_loc;
            accum_wtime_sum[tid] = wtime_sum_loc;

         }

         //printf("%d %e %e\n", tid, read_wtime_sum[tid], write_wtime_sum[tid]);

         //PrintDummy(dummy);
      }

      /* Compute means of all measurements */
      uint64_t read_cycles_mean = 0;
      for (int t = 0; t < num_threads; t++){
         read_cycles_mean += read_cycles_sum[t];
      }
      read_cycles_mean /= (uint64_t)num_threads_read;
      double read_wtime_mean = 0;
      for (int t = 0; t < num_threads; t++){
         read_wtime_mean += read_wtime_sum[t];
      }
      read_wtime_mean /= (double)num_threads_read;

      uint64_t write_cycles_mean = 0;
      for (int t = 0; t < num_threads; t++){
         write_cycles_mean += write_cycles_sum[t];
      }
      write_cycles_mean /= (uint64_t)num_threads_write;
      double write_wtime_mean = 0;
      for (int t = 0; t < num_threads; t++){
         write_wtime_mean += write_wtime_sum[t];
      }
      write_wtime_mean /= (double)num_threads_write;

      uint64_t accum_cycles_mean = 0;
      for (int t = 0; t < num_threads; t++){
         accum_cycles_mean += accum_cycles_sum[t];
      }
      accum_cycles_mean /= (uint64_t)num_threads;
      double accum_wtime_mean = 0;
      for (int t = 0; t < num_threads; t++){
         accum_wtime_mean += accum_wtime_sum[t];
      }
      accum_wtime_mean /= (double)num_threads;

      /* Print output */
      if (verbose_output){
         printf("Read wall-clock time  = %e\n"
                "Read num cycles = %" PRIu64 "\n"
                "Write wall-clock time  = %e\n"
                "Write num cycles = %" PRIu64 "\n"
                "Accum wall-clock time  = %e\n"
                "Accum num cycles = %" PRIu64 "\n"
                "Acummulated value = %e\n", 
                read_wtime_mean, read_cycles_mean,
                write_wtime_mean, write_cycles_mean,
                accum_wtime_mean, accum_cycles_mean,
                a);
      }
      else {
         printf("%e %" PRIu64 " "
                "%e %" PRIu64 " "
                "%e %" PRIu64 " "
                "%e\n",
                read_wtime_mean, read_cycles_mean,
                write_wtime_mean, write_cycles_mean,
                accum_wtime_mean, accum_cycles_mean,
                a);
      }
   }

   free(read_cycles_sum);
   free(read_wtime_sum);
   free(write_cycles_sum);
   free(write_wtime_sum);
   free(accum_cycles_sum);
   free(accum_wtime_sum);
   free(a);

   return 0;
}
