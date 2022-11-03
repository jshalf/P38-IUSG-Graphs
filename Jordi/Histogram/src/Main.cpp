#include "Main.hpp"
#include "Histogram.hpp"
#include "Matrix.hpp"
#include "Misc.hpp"

int main (int argc, char *argv[])
{
   /* matrix defaults */
   SparseMatrixInput spm_input;
   spm_input.mat_type = MatrixType::lower;
   spm_input.file_type = FileType::bin;
   spm_input.num_rows = 100;
   spm_input.num_cols = 100;
   spm_input.grid_len.resize(1);
   spm_input.grid_len[0] = 10;
   spm_input.max_row_nnz = 5;
   spm_input.store_type = SparseMatrixStorageType::CSR;
   spm_input.store_diag_in_vec = true;

   HistogramData hist;
   hist.input.num_threads = 1;
   hist.input.atomic_flag = 1;
   hist.input.MsgQ_flag = 0;
   hist.input.comp_wtime_flag = 0;
   hist.input.MsgQ_wtime_flag = 0;
   hist.input.MsgQ_cycles_flag = 0;
   hist.input.comp_noop_flag = 0;
   hist.input.MsgQ_noop_flag = 0;
   hist.input.reduce_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   int problem_type = PROBLEM_5PT_POISSON;
   char mat_file_str[128];

   int print_usage = 0;
   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem*/
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-algo") == 0){ /* method */
         arg_index++;
         if (strcmp(argv[arg_index], "async") == 0){ /* async method that uses atomics or msgQs */
            hist.input.reduce_flag = 0;
         }
         else if (strcmp(argv[arg_index], "reduce") == 0){ /* each thread uses a local vector which are summed at the end */
            hist.input.reduce_flag = 1;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* ``size'' of matrix. n*n rows for Laplace, n rows otherwise */
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         hist.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use no-atomics implementations */
         hist.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of TriSolve runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues */
         hist.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_wtime") == 0){
         hist.input.MsgQ_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_wtime") == 0){
         hist.input.comp_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_cycles") == 0){
         hist.input.MsgQ_cycles_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_noop") == 0){
         hist.input.MsgQ_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_noop") == 0){
         hist.input.comp_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-problem <problem_name>:  test problem.\n");
      printf("      5pt:                five-point centered difference discretization of the Poisson equation.\n");
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-method <method_name>:    method for computing y=A^Tx.\n");
      printf("      atomic:             method that uses atomics.\n");
      printf("      reduce:             baseline method that does not require atomics.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("\n");
      return 0;
   }
   
   omp_set_num_threads(hist.input.num_threads);

   int csc_flag = 0, coo_flag = 0;
   int include_diag = 1;

   /* set up problem */
   SparseMatrix A(spm_input);

   if (problem_type == PROBLEM_FILE){
      A.ConstructMatrixFromFile();
   }
   else if (problem_type == PROBLEM_5PT_POISSON){
      A.ConstructLaplace2D5pt();
   }
   else {
      srand(0);
      A.ConstructRandomMatrix();
   }

   vector<int> index = A.GetColIndices();
   int n_index = A.GetNNZ();
   int n_tally = A.GetNumRows();

   int *Tally_exact = (int *)calloc(n_tally, sizeof(int));
   int *Tally = (int *)calloc(n_tally, sizeof(int));
   double *e = (double *)calloc(n_tally, sizeof(double));
   double *y = (double *)calloc(n_tally, sizeof(double));

   if (hist.input.expand_flag == 1){
      hist.Tally_expand = (int *)calloc(hist.input.num_threads * n_tally, sizeof(int));
   }

   hist.output.solve_wtime_vec = (double *)calloc(hist.input.num_threads, sizeof(double));
   hist.output.num_qGets_vec = (int *)calloc(hist.input.num_threads, sizeof(int));
   hist.output.num_qPuts_vec = (int *)calloc(hist.input.num_threads, sizeof(int));
  
   /* output stats for msgQ stuff */ 
   if (hist.input.MsgQ_flag == 1){
      //if (hist.input.MsgQ_wtime_flag == 1){
         hist.output.MsgQ_put_wtime_vec = (double *)calloc(hist.input.num_threads, sizeof(double));
         hist.output.MsgQ_get_wtime_vec = (double *)calloc(hist.input.num_threads, sizeof(double));
      //}
      //else if (hist.input.MsgQ_cycles_flag == 1){
         hist.output.MsgQ_put_cycles_vec = (uint64_t *)calloc(hist.input.num_threads, sizeof(uint64_t));
         hist.output.MsgQ_get_cycles_vec = (uint64_t *)calloc(hist.input.num_threads, sizeof(uint64_t));
      //}
      //else if (hist.input.comp_wtime_flag == 1){
         hist.output.comp_wtime_vec = (double *)calloc(hist.input.num_threads, sizeof(double));
      //}
   }

   srand(0);

   /* do the sequential hist first,
      we'll use this to check that the parallel hist is correct */
   Histogram_Seq(&hist, index, n_index, Tally_exact, n_tally);
   //#pragma omp parallel
   //{
   //   double dummy = 0.0;
   //   #pragma omp for schedule(static, 1) nowait
   //   for (int i = 0; i < n_index; i++){
   //      dummy += index[i];
   //   }
   //   PrintDummy(dummy);
   //}

   for (int run = 1; run <= num_runs; run++){
      for (int i = 0; i < n_tally; i++){
         Tally[i] = 0;
      }

      /* parallel Histogram */
      Histogram_Par(&hist, index, n_index, Tally, n_tally);

      /* compute error */
      for (int i = 0; i < n_tally; i++){
         y[i] = (double)Tally_exact[i];
         e[i] = (double)(Tally_exact[i] - Tally[i]);
         //printf("%d %d\n", Tally_exact[i], Tally[i]);
      }

      double error_norm = sqrt(InnerProd(e, e, n_tally))/sqrt(InnerProd(y, y, n_tally));

      double solve_wtime_sum = accumulate(hist.output.solve_wtime_vec, hist.output.solve_wtime_vec+hist.input.num_threads, (double)0.0);
      double solve_wtime_mean = solve_wtime_sum / (double)hist.input.num_threads;

      /* compute msgQ output stats */
      double MsgQ_put_wtime_sum = 0.0, MsgQ_put_wtime_mean = 0.0;
      double MsgQ_get_wtime_sum = 0.0, MsgQ_get_wtime_mean = 0.0;
      uint64_t MsgQ_put_cycles_sum = 0;
      uint64_t MsgQ_get_cycles_sum = 0;
      double MsgQ_put_cycles_mean = 0.0;
      double MsgQ_get_cycles_mean = 0.0;
      double comp_wtime_sum = 0.0, comp_wtime_mean = 0.0;
      int num_qGets_sum = 0, num_qPuts_sum = 0;
      if (hist.input.MsgQ_flag == 1){
         //if (hist.input.MsgQ_wtime_flag == 1){
            MsgQ_put_wtime_sum = accumulate(hist.output.MsgQ_put_wtime_vec, hist.output.MsgQ_put_wtime_vec+hist.input.num_threads, (double)0.0);
            MsgQ_put_wtime_mean = MsgQ_put_wtime_sum / (double)hist.input.num_threads;

            MsgQ_get_wtime_sum = accumulate(hist.output.MsgQ_get_wtime_vec, hist.output.MsgQ_get_wtime_vec+hist.input.num_threads, (double)0.0);
            MsgQ_get_wtime_mean = MsgQ_get_wtime_sum / (double)hist.input.num_threads;
         //}
         //else if (hist.input.MsgQ_cycles_flag == 1){
            MsgQ_put_cycles_sum = accumulate(hist.output.MsgQ_put_cycles_vec, hist.output.MsgQ_put_cycles_vec+hist.input.num_threads, (uint64_t)0);
            MsgQ_put_cycles_mean = MsgQ_put_cycles_sum / (double)hist.input.num_threads;
            
            MsgQ_get_cycles_sum = accumulate(hist.output.MsgQ_get_cycles_vec, hist.output.MsgQ_get_cycles_vec+hist.input.num_threads, (uint64_t)0);
            MsgQ_get_cycles_mean = MsgQ_get_cycles_sum / (double)hist.input.num_threads;
         //}
         //else if (hist.input.comp_wtime_flag == 1){
            comp_wtime_sum = accumulate(hist.output.comp_wtime_vec, hist.output.comp_wtime_vec+hist.input.num_threads, 0.0);
            comp_wtime_mean = comp_wtime_sum / (double)hist.input.num_threads;
         //}
         
         num_qGets_sum = accumulate(hist.output.num_qGets_vec, hist.output.num_qGets_vec+hist.input.num_threads, (int)0);
         num_qPuts_sum = accumulate(hist.output.num_qPuts_vec, hist.output.num_qPuts_vec+hist.input.num_threads, (int)0);
      }

      /* print stats */
      if (verbose_output){
         printf("Error L2-norm = %e\n"
                "Wall-clock time = %e\n"
                "Mean comp wtime = %e\n"
                "MsgQ put wtime = %e\n"
                "MsgQ get wtime = %e\n"
                "MsgQ put cycles = %" PRIu64 "\n"
                "MsgQ get cycles = %" PRIu64 "\n"
                "Num qPuts = %d\n"
                "Num qGets = %d\n",
                error_norm,
                solve_wtime_mean,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
      else {
         printf("%e %e %e %e %e %" PRIu64 " %" PRIu64 " %d %d\n",
                error_norm,
                solve_wtime_mean,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
   }

   free(Tally_exact);
   free(Tally);
   free(y);
   free(e);

   if (hist.input.expand_flag == 1){
      free(hist.Tally_expand);
   }

   free(hist.output.solve_wtime_vec);
   free(hist.output.num_qGets_vec);
   free(hist.output.num_qPuts_vec);

   if (hist.input.MsgQ_flag == 1){
      //if (hist.input.MsgQ_wtime_flag == 1){
         free(hist.output.MsgQ_put_wtime_vec);
         free(hist.output.MsgQ_get_wtime_vec);
      //}
      //else if (hist.input.MsgQ_cycles_flag == 1){
         free(hist.output.MsgQ_put_cycles_vec);
         free(hist.output.MsgQ_get_cycles_vec);
      //}
      //else if (hist.input.comp_wtime_flag == 1){
         free(hist.output.comp_wtime_vec);
      //}
   }

   return 0;
}
