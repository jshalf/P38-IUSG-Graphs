#include "Main.hpp"
#include "Misc.hpp"
#include "MsgQ.hpp"

using namespace std;

typedef struct {
   int src_row;
   int dst_row;
   double data;
} PingPongMessage;

int main (int argc, char *argv[])
{
#ifdef USE_DEVA
   int deva_flag = 1;
#else
   int deva_flag = 0;
#endif
   int num_threads = 1;
   int verbose_output = 0;
   int num_runs = 1;
   int num_iters = 100;
   int num_rounds = 1;
   int stride = 1;
   int concurrent_comm = 0;

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-num_iters") == 0){
         arg_index++;
         num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_rounds") == 0){
         arg_index++;
         num_rounds = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-comm_scheme") == 0){
         arg_index++;
         if (strcmp(argv[arg_index], "two") == 0){
         }
         else if (strcmp(argv[arg_index], "all") == 0){
         }
         else if (strcmp(argv[arg_index], "rand") == 0){
         }
      }
      else if (strcmp(argv[arg_index], "-contention") == 0){
         concurrent_comm = 1;
      }
      else if (strcmp(argv[arg_index], "-stride") == 0){
         arg_index++;
         stride = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-num_threads <int value>:       Number of threads (OpenMP or devastator).\n");
      printf("-num_runs <int value>:          Number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                      Verbose output.\n");
      printf("-num_iters <int value>:         Number of ping-pong experiments\n");
      printf("-num_rounds <int value>:        Number of times a message is sent back and forth\n");
      printf("\n");
      return 0;
   }

   srand(time(NULL));
  
#ifdef USE_DEVA 
#else
   omp_set_num_threads(num_threads);
#endif


#ifdef USE_DEVA
   MessageQueue<PingPongMessage> *Q;
   Q = new MessageQueue<PingPongMessage>(stride * num_threads);
#else
   vector<MessageQueue<PingPongMessage>*> Q;
   Q.resize(stride * num_threads);
   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      Q[stride * tid] = new MessageQueue<PingPongMessage>(stride * num_threads);
   }
#endif

   vector<vector<vector<uint64_t>>> put_rounds_cycles_sum(
                                       num_rounds,
                                       vector<vector<uint64_t>>(num_threads,
                                                                vector<uint64_t>(num_threads,
                                                                                 (uint64_t)0))
                                                               );
   vector<vector<vector<uint64_t>>> get_rounds_cycles_sum(
                                       num_rounds,
                                       vector<vector<uint64_t>>(num_threads,
                                                                vector<uint64_t>(num_threads,
                                                                                 (uint64_t)0))
                                                               );
   vector<vector<vector<uint64_t>>> inflight_rounds_cycles_sum(
                                       num_rounds,
                                       vector<vector<uint64_t>>(num_threads,
                                                                vector<uint64_t>(num_threads,
                                                                                 (uint64_t)0))
                                                               );

   vector<vector<vector<int>>> put_rounds_measure_count(
                                       num_rounds,
                                       vector<vector<int>>(num_threads,
                                                                vector<int>(num_threads,
                                                                                 (int)0))
                                                       );
   vector<vector<vector<int>>> get_rounds_measure_count(
                                       num_rounds,
                                       vector<vector<int>>(num_threads,
                                                                vector<int>(num_threads,
                                                                                 (int)0))
                                                       );
   vector<vector<vector<int>>> inflight_rounds_measure_count(
                                       num_rounds,
                                       vector<vector<int>>(num_threads,
                                                                vector<int>(num_threads,
                                                                                 (int)0))
                                                            );
 
   for (int run = 1; run <= num_runs; run++){
      for (int i = 0; i < num_threads; i++){
         for (int j = 0; j < num_threads; j++){
            for (int round = 0; round < num_rounds; round++){
               put_rounds_cycles_sum[round][i][j] = (uint64_t)0;
               get_rounds_cycles_sum[round][i][j] = (uint64_t)0;
               inflight_rounds_cycles_sum[round][i][j] = (uint64_t)0;
               put_rounds_measure_count[round][i][j] = 0;
               get_rounds_measure_count[round][i][j] = 0;
               inflight_rounds_measure_count[round][i][j] = 0;
            }
         }
      }

      double b_glob = RandDouble(-1.0, 1.0);
      for (int iter = 0; iter < num_iters; iter++){
         for (int put_tid = 0; put_tid < num_threads; put_tid++){
            for (int get_tid = 0; get_tid < num_threads; get_tid++){
               if (put_tid != get_tid){
#ifdef USE_DEVA
                  deva::run( [&Q,
                              &b_glob,
                              num_threads,
                              stride,
                              iter,
                              num_iters,
                              num_rounds,
                              put_tid,
                              get_tid,
                              concurrent_comm,
                              &put_rounds_cycles_sum,
                              &get_rounds_cycles_sum,
                              &inflight_rounds_cycles_sum,
                              &put_rounds_measure_count,
                              &get_rounds_measure_count,
                              &inflight_rounds_measure_count] () {
                     int tid = deva::rank_me();
#else
                  #pragma omp parallel
                  {
                     int tid = omp_get_thread_num();
#endif
                     double b = RandDouble(-1.0, 1.0);

                     vector<vector<uint64_t>> put_rounds_cycles_sum_loc(num_rounds,
                                                                        vector<uint64_t>(num_threads, (uint64_t)0));
                     vector<vector<uint64_t>> get_rounds_cycles_sum_loc(num_rounds,
                                                                        vector<uint64_t>(num_threads, (uint64_t)0));
                     vector<vector<uint64_t>> inflight_rounds_cycles_sum_loc(num_rounds,
                                                                             vector<uint64_t>(num_threads, (uint64_t)0));

                     vector<vector<int>> put_rounds_measure_count_loc(num_rounds, vector<int>(num_threads, (int)0));
                     vector<vector<int>> get_rounds_measure_count_loc(num_rounds, vector<int>(num_threads, (int)0));
                     vector<vector<int>> inflight_rounds_measure_count_loc(num_rounds, vector<int>(num_threads, (int)0));

                     int put_round = 0;
                     int get_round = 0;
                     int temp_put_tid, temp_get_tid, tot_rounds;
                     if (concurrent_comm == 1){
                        temp_put_tid = tid;
                        temp_get_tid = (put_tid == tid ? get_tid : put_tid);
                        tot_rounds = num_rounds;
                     }
                     else {
                        temp_put_tid = put_tid;
                        temp_get_tid = get_tid;
                        tot_rounds = 2 * num_rounds;
                     }

                     for (int round = 1; round <= tot_rounds; round++){
                        //if (tid == 0) printf("%d %d %d %d\n", iter, round, _tid, temp_get_tid);
                        if (tid == temp_put_tid){
                           double a = RandDouble(-1.0, 1.0);
                           uint64_t cycles_start = rdtsc();
                           PingPongMessage msg;
                           msg.src_row = tid;
                           msg.dst_row = temp_get_tid;
                           msg.data = a;
#ifdef USE_DEVA
                           Q->qPut(stride * temp_get_tid, msg);
                           //deva::progress();
#else
                           Q[stride * temp_get_tid]->qPut(stride * tid, msg);
#endif
                           uint64_t cycles_stop = rdtsc();
                           uint64_t cycles_elapsed = cycles_stop - cycles_start;

                           put_rounds_cycles_sum_loc[put_round][temp_get_tid] += cycles_elapsed;
                           put_rounds_measure_count_loc[put_round][temp_get_tid]++;

                           put_round++;
                        }
                        if (concurrent_comm == 1){
                           temp_put_tid = temp_get_tid;
                           temp_get_tid = tid;
                        }
                        if (tid == temp_get_tid){
                           uint64_t cycles_start = rdtsc();
                           do {
                              PingPongMessage msg;

                              uint64_t inner_cycles_start = rdtsc();
#ifdef USE_DEVA
                              deva::progress();
                              int get_flag = Q->qGet(stride * tid, &msg);
#else
                              int get_flag = Q[stride * tid]->qGet(stride * temp_put_tid, &msg);
#endif
                              if (get_flag){
                                 uint64_t cycles_stop = rdtsc();
                                 uint64_t inner_cycles_elapsed = cycles_stop - inner_cycles_start;
                                 uint64_t cycles_elapsed = cycles_stop - cycles_start;

                                 get_rounds_cycles_sum_loc[get_round][temp_put_tid] += inner_cycles_elapsed;
                                 get_rounds_measure_count_loc[get_round][temp_put_tid]++;

                                 inflight_rounds_cycles_sum_loc[get_round][temp_put_tid] +=
                                    (cycles_elapsed - inner_cycles_elapsed);
                                 inflight_rounds_measure_count_loc[get_round][temp_put_tid]++;
                                 
                                 get_round++;                                 

                                 //b += msg.data + (double)msg.src_row + (double)msg.dst_row;
                                 break;
                              } 
                           } while (1);
                        }
                        int temp = temp_put_tid;
                        temp_put_tid = temp_get_tid;
                        temp_get_tid = temp;
                     }

#ifdef USE_DEVA
#else
                     //#pragma omp barrier
                     //#pragma omp atomic
#endif
                     //b_glob += b;


                     if (iter == num_iters-1){
                        //for (int t = 0; t < num_threads; t++){
                        //   put_cycles_sum[tid][t] += put_cycles_sum_loc[t];
                        //   get_cycles_sum[tid][t] += get_cycles_sum_loc[t];
                        //   inflight_cycles_sum[tid][t] += inflight_cycles_sum_loc[t];
                        //   
                        //   put_measure_count[tid][t] += put_measure_count_loc[t];
                        //   get_measure_count[tid][t] += get_measure_count_loc[t];
                        //   inflight_measure_count[tid][t] += inflight_measure_count_loc[t];
                        //}
                        for (int round = 0; round < num_rounds; round++){
                           for (int t = 0; t < num_threads; t++){
                              put_rounds_cycles_sum[round][tid][t] += put_rounds_cycles_sum_loc[round][t];
                              get_rounds_cycles_sum[round][tid][t] += get_rounds_cycles_sum_loc[round][t];
                              inflight_rounds_cycles_sum[round][tid][t] += inflight_rounds_cycles_sum_loc[round][t];

                              put_rounds_measure_count[round][tid][t] += put_rounds_measure_count_loc[round][t];
                              get_rounds_measure_count[round][tid][t] += get_rounds_measure_count_loc[round][t];
                              inflight_rounds_measure_count[round][tid][t] += inflight_rounds_measure_count_loc[round][t];
                           }
                        }
                     }
#ifdef USE_DEVA
                  });
#else
                  }
#endif
               }
            }
         }
      }

      
      //for (int i = 0; i < num_threads; i++){
      //   for (int j = 0; j < num_threads; j++){
      //      double put_cycles_mean = (double)put_cycles_sum[i][j];
      //      put_cycles_mean /= max(1.0, (double)put_measure_count[i][j]);

      //      double get_cycles_mean = (double)get_cycles_sum[i][j];
      //      get_cycles_mean /= max(1.0, (double)get_measure_count[i][j]);

      //      double inflight_cycles_mean = (double)inflight_cycles_sum[i][j];
      //      inflight_cycles_mean /= max(1.0, (double)inflight_measure_count[i][j]);

      //      printf("%d %d "
      //             "%f "
      //             "%f "
      //             "%f\n",
      //             i, j, 
      //             put_cycles_mean,
      //             get_cycles_mean,
      //             inflight_cycles_mean);
      //   }
      //}

      for (int round = 0; round < num_rounds; round++){
         for (int i = 0; i < num_threads; i++){
            for (int j = 0; j < num_threads; j++){
               double put_cycles_mean = (double)put_rounds_cycles_sum[round][i][j];
               put_cycles_mean /= max(1.0, (double)put_rounds_measure_count[round][i][j]);

               double get_cycles_mean = (double)get_rounds_cycles_sum[round][i][j];
               get_cycles_mean /= max(1.0, (double)get_rounds_measure_count[round][i][j]);

               double inflight_cycles_mean = (double)inflight_rounds_cycles_sum[round][i][j];
               inflight_cycles_mean /= max(1.0, (double)inflight_rounds_measure_count[round][i][j]);

               //if (i == 31 && j == 0)
               printf("%d %d %d "
                      "%f "
                      "%f "
                      "%f\n",
                      round, i, j,
                      put_cycles_mean,
                      get_cycles_mean,
                      inflight_cycles_mean);
            }
         }
      }

      char dummy_filename[128];
      sprintf(dummy_filename, "dummy.txt");
      FILE *file = fopen(dummy_filename, "w");
      fprintf(file, "%e\n", b_glob);
      fclose(file);

      //uint64_t get_cycles_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   get_cycles_mean += get_cycles_sum[t];
      //}
      //get_cycles_mean /= (uint64_t)num_threads_get;
      //double get_wtime_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   get_wtime_mean += get_wtime_sum[t];
      //}
      //get_wtime_mean /= (double)num_threads_get;

      //uint64_t put_cycles_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   put_cycles_mean += put_cycles_sum[t];
      //}
      //put_cycles_mean /= (uint64_t)num_threads_put;
      //double put_wtime_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   put_wtime_mean += put_wtime_sum[t];
      //}
      //put_wtime_mean /= (double)num_threads_put;

      //uint64_t accum_cycles_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   accum_cycles_mean += accum_cycles_sum[t];
      //}
      //accum_cycles_mean /= (uint64_t)num_threads;
      //double accum_wtime_mean = 0;
      //for (int t = 0; t < num_threads; t++){
      //   accum_wtime_mean += accum_wtime_sum[t];
      //}
      //accum_wtime_mean /= (double)num_threads;

      //if (verbose_output){
      //   printf("Read wall-clock time  = %e\n"
      //          "Read num cycles = %" PRIu64 "\n"
      //          "Write wall-clock time  = %e\n"
      //          "Write num cycles = %" PRIu64 "\n"
      //          "Accum wall-clock time  = %e\n"
      //          "Accum num cycles = %" PRIu64 "\n"
      //          "Acummulated value = %e\n", 
      //          get_wtime_mean, get_cycles_mean,
      //          put_wtime_mean, put_cycles_mean,
      //          accum_wtime_mean, accum_cycles_mean,
      //          a);
      //}
      //else {
      //   printf("%e %" PRIu64 " "
      //          "%e %" PRIu64 " "
      //          "%e %" PRIu64 " "
      //          "%e\n",
      //          get_wtime_mean, get_cycles_mean,
      //          put_wtime_mean, put_cycles_mean,
      //          accum_wtime_mean, accum_cycles_mean,
      //          a);
      //}
   }

   return 0;
}
