#include "Main.hpp"
#include "Matrix.hpp"
#include "Parallel.hpp"
#include "Misc.hpp"
#include "MsgQ.hpp"
#include "LevelSets.hpp"
#include "TriSolve.hpp"

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

   /* set defaults */
   TriSolverSolveOrder tri_solve_order = TriSolverSolveOrder::natural;
   CommunicationType comm_type = CommunicationType::atomic;

   TriSolveData ts;
   ts.input.num_threads = 1;
   ts.input.atomic_flag = 1;
   ts.input.async_flag = 1;
   ts.input.omp_for_flag = 1;
   ts.input.MsgQ_flag = 0;
   ts.input.block_size = 1;
   ts.input.fine_grained_flag = 0;
   ts.input.comp_wtime_flag = 0;
   ts.input.MsgQ_wtime_flag = 0;
   ts.input.MsgQ_cycles_flag = 0;
   ts.input.comp_noop_flag = 0;
   ts.input.MsgQ_noop_flag = 0;
   ts.input.solver_type = TRISOLVE_ASYNC;
   ts.input.setup_type = LEVEL_SETS_SEQ_SETUP;
   int verbose_output = 0;
   int num_runs = 1;
   int problem_type = PROBLEM_RANDOM;
   double start;

   /* command line arguments */
   int arg_index = 0;
   int print_usage = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem*/
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "rand") == 0){ /* random matrix */
            problem_type = PROBLEM_RANDOM;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(spm_input.file_name, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-solver") == 0){ /* solver type */
         arg_index++;
         if (strcmp(argv[arg_index], "async") == 0){ /* asynchronous */
            ts.input.solver_type = TRISOLVE_ASYNC;
            tri_solve_order = TriSolverSolveOrder::natural;
         }
         else if (strcmp(argv[arg_index], "lev_sched") == 0){ /* classical level-scheduled */
            ts.input.solver_type = TRISOLVE_LEVEL_SCHEDULED;
            tri_solve_order = TriSolverSolveOrder::level_sched;
         }
         else if (strcmp(argv[arg_index], "async_lev_sched") == 0){ /* asynchronous using level-set info */
            ts.input.solver_type = TRISOLVE_ASYNC_LEVEL_SCHEDULED;
            tri_solve_order = TriSolverSolveOrder::level_sched;
         }
      }
      else if (strcmp(argv[arg_index], "-setup") == 0){ /* solver type */
         arg_index++;
         if (strcmp(argv[arg_index], "async") == 0){ /* asynchronous */
            ts.input.setup_type = LEVEL_SETS_ASYNC_SETUP;
         }
         else if (strcmp(argv[arg_index], "seq") == 0){ /* classical level-scheduled */
            ts.input.setup_type = LEVEL_SETS_SEQ_SETUP;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* ``size'' of matrix. n*n rows for Laplace, n rows otherwise. */
         arg_index++;
         spm_input.num_rows = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         ts.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use no-atomics implementations */
         ts.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of TriSolve runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues in async solvers */
         ts.input.MsgQ_flag = 1;
         comm_type = CommunicationType::MsgQ;
      }
      else if (strcmp(argv[arg_index], "-mxr_nnz") == 0){ /* max number of non-zeros per row (only used for generating random amtrices) */
         arg_index++;
         spm_input.max_row_nnz = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-sync") == 0){ /* make the async solver synchronous */
         ts.input.async_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-sp_store_type") == 0){ /* matrix storage type */
         arg_index++;
         if (strcmp(argv[arg_index], "csr") == 0){ /* compressed sparse row */
            spm_input.store_type = SparseMatrixStorageType::CSR;
         }
         else if (strcmp(argv[arg_index], "csc") == 0){ /* compressed sparse column */
            spm_input.store_type = SparseMatrixStorageType::CSC;
         }
      }
      else if (strcmp(argv[arg_index], "-fine_grained") == 0){ /* use fine grained solver */
         ts.input.fine_grained_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_wtime") == 0){
         ts.input.MsgQ_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_wtime") == 0){
         ts.input.comp_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_cycles") == 0){
         ts.input.MsgQ_cycles_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_noop") == 0){
         ts.input.MsgQ_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_noop") == 0){
         ts.input.comp_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-block_size") == 0){
         arg_index++;
         ts.input.block_size = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-problem <string>:  test problem.\n");
      printf("      5pt:                five-point centered difference discretization of the Poisson equation.\n");
      printf("      rand:               random sparse matrix.\n");
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-solver <string>:    method for solving LUx=b.\n");
      printf("      lev_sched:          level scheduled method.\n");
      printf("      async:              asynchronous method.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("                          for the random matrix, this is the number of rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measurements and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("-mxr_nnz <int value>:     maximum number of non-zero values per row for the random matrix.\n");
      printf("-sp_store_type <string>:  type of sparse matrix storage.\n");
      printf("      csr:          compressed sparse row.\n");
      printf("      csc:          compressed sparse column.\n");
      printf("-fine_grained:            fine-grained version of asynchronous solver.\n");
      printf("\n");
      return 0;
   }

   SparseMatrix L(spm_input);

   if (problem_type == PROBLEM_FILE){
      L.ConstructMatrixFromFile();
   }
   else { 
      srand(0);
      L.ConstructRandomMatrix();
   }

   //char L_outfile[128];
   //if (csc_flag == 1){
   //   sprintf(L_outfile, "./matlab/L_csc.txt");
   //}
   //else {
   //   sprintf(L_outfile, "./matlab/L_csr.txt");
   //}
   //PrintMatrix(L, L_outfile, 1, csc_flag);

   int num_rows = L.GetNumRows();
   
   ts.output.setup_wtime = 0.0;
   /* set up stuff for serial solver */
   int *L_perm = (int *)calloc(num_rows, sizeof(int));
   //int *U_perm = (int *)calloc(num_rows, sizeof(int));
   for (int i = 0; i < num_rows; i++){
      L_perm[i] = i;
      //U_perm[i] = num_rows - (i+1);
   }

   std::vector<double> x(num_rows), x_exact(num_rows), e_x(num_rows), b(num_rows);

   srand(0);
   for (int i = 0; i < num_rows; i++){
      b[i] = 1.0;//RandDouble(-1.0, 1.0);
   }

   ts.output.row_output = (RowOutputData *)malloc(num_rows * sizeof(RowOutputData));
   for (int i = 0; i < num_rows; i++){
      ts.output.row_output[i].wtime = 0.0;
      ts.output.row_output[i].num_spins = 0;
      ts.output.row_output[i].num_atomics = 0;
   }

   ts.output.setup_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
   ts.output.solve_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
   ts.output.num_relax = (int *)calloc(ts.input.num_threads, sizeof(int));
   ts.output.num_iters = (int *)calloc(ts.input.num_threads, sizeof(int));

   if (ts.input.MsgQ_flag == 1){
      //if (ts.input.MsgQ_wtime_flag == 1){
         ts.output.MsgQ_put_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
         ts.output.MsgQ_get_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
      //}
      //else if (ts.input.MsgQ_cycles_flag == 1){
         ts.output.MsgQ_put_cycles_vec = (uint64_t *)calloc(ts.input.num_threads, sizeof(uint64_t));
         ts.output.MsgQ_get_cycles_vec = (uint64_t *)calloc(ts.input.num_threads, sizeof(uint64_t));
      //}
      //else if (ts.input.comp_wtime_flag == 1){
         ts.output.comp_wtime_vec = (double *)calloc(ts.input.num_threads, sizeof(double));
      //}
      
      ts.output.num_qGets_vec = (int *)calloc(ts.input.num_threads, sizeof(int));
      ts.output.num_qPuts_vec = (int *)calloc(ts.input.num_threads, sizeof(int));
   }

   LevelSchedTriSolver seq_tri_solver(L, 1);
   LevelSchedTriSolver *lev_sched_tri_solver;
   AsyncTriSolver *async_tri_solver;

#if USE_DEVA
#elif  USE_STDTHREADS
#else
   omp_set_num_threads(1);
#endif
   /* serial solver first */
   double seq_start = omp_get_wtime();
   seq_tri_solver.Setup(ts.output);
   seq_tri_solver.Solve(b, x_exact, ts.output);
   double seq_wtime = omp_get_wtime() - seq_start;
   double init_setup_time = 0;


   int lvl_n_min = 0, lvl_n_max = 0, num_lvls = 0;
   double lvl_n_mean = 0.0;
#if USE_DEVA
   ts.input.num_threads = proc_rank_n;
#elif  USE_STDTHREADS
#else
   omp_set_num_threads(ts.input.num_threads);
#endif
   if (ts.input.solver_type == TRISOLVE_LEVEL_SCHEDULED){
      lev_sched_tri_solver = new LevelSchedTriSolver(L, ts.input.num_threads);
      lev_sched_tri_solver->Setup(ts.output);

      //lvl_n_min = *min_element(ts.L_lvl_set.level_size.begin(), ts.L_lvl_set.level_size.end());
      //lvl_n_max = *max_element(ts.L_lvl_set.level_size.begin(), ts.L_lvl_set.level_size.end());
      //lvl_n_mean = (double)accumulate(ts.L_lvl_set.level_size.begin(), ts.L_lvl_set.level_size.end(), (int)0) / ts.L_lvl_set.level_size.size(); 
      //num_lvls = ts.L_lvl_set.level_size.size();
   }
   else {
      async_tri_solver = new AsyncTriSolver(L, ts.input.num_threads);
      async_tri_solver->Setup(ts.output, tri_solve_order, comm_type); 
   }

   for (int run = 1; run <= num_runs; run++){
      ts.output.setup_wtime = init_setup_time; 
      if (ts.input.solver_type == TRISOLVE_LEVEL_SCHEDULED){
         lev_sched_tri_solver->InitSolveData();
      }
      else {
         async_tri_solver->InitSolveData();
      }

      double start = omp_get_wtime();
      if (ts.input.solver_type == TRISOLVE_LEVEL_SCHEDULED){
         lev_sched_tri_solver->Solve(b, x, ts.output);
      }
      else {
         async_tri_solver->Solve(b, x, ts.output);
      }
      double overall_solve_wtime = omp_get_wtime() - start;

      double setup_wtime_sum = SumDouble(ts.output.setup_wtime_vec, ts.input.num_threads);
      ts.output.setup_wtime += setup_wtime_sum / (double)ts.input.num_threads;

      double solve_wtime_sum = SumDouble(ts.output.solve_wtime_vec, ts.input.num_threads);
      ts.output.solve_wtime = solve_wtime_sum / (double)ts.input.num_threads;

      /* compute the error between the serial and parallel solvers */
      for (int i = 0; i < num_rows; i++){
         e_x[i] = x_exact[i] - x[i];
         //printf("%e %e\n", x_exact[i], x[i]);
      }
      double error_x = sqrt(InnerProd(e_x.data(), e_x.data(), num_rows)) /
                       sqrt(InnerProd(x_exact.data(), x_exact.data(), num_rows));
      int num_relax_sum = SumInt(ts.output.num_relax, ts.input.num_threads);
      int num_iters_sum = SumInt(ts.output.num_iters, ts.input.num_threads);

      double MsgQ_put_wtime_sum = 0.0, MsgQ_put_wtime_mean = 0.0;
      double MsgQ_get_wtime_sum = 0.0, MsgQ_get_wtime_mean = 0.0;
      uint64_t MsgQ_put_cycles_sum = 0;
      uint64_t MsgQ_get_cycles_sum = 0;
      double MsgQ_put_cycles_mean = 0.0;
      double MsgQ_get_cycles_mean = 0.0;
      double comp_wtime_sum = 0.0, comp_wtime_mean = 0.0;
      int num_qGets_sum = 0, num_qPuts_sum = 0;
      if (ts.input.MsgQ_flag == 1){
         //if (ts.input.MsgQ_wtime_flag == 1){
            MsgQ_put_wtime_sum = accumulate(ts.output.MsgQ_put_wtime_vec, ts.output.MsgQ_put_wtime_vec+ts.input.num_threads, (double)0.0);
            MsgQ_put_wtime_mean = MsgQ_put_wtime_sum / (double)ts.input.num_threads;

            MsgQ_get_wtime_sum = accumulate(ts.output.MsgQ_get_wtime_vec, ts.output.MsgQ_get_wtime_vec+ts.input.num_threads, (double)0.0);
            MsgQ_get_wtime_mean = MsgQ_get_wtime_sum / (double)ts.input.num_threads;
         //}
         //else if (ts.input.MsgQ_cycles_flag == 1){
            MsgQ_put_cycles_sum = accumulate(ts.output.MsgQ_put_cycles_vec, ts.output.MsgQ_put_cycles_vec+ts.input.num_threads, (uint64_t)0);
            MsgQ_put_cycles_mean = MsgQ_put_cycles_sum / (double)ts.input.num_threads;

            MsgQ_get_cycles_sum = accumulate(ts.output.MsgQ_get_cycles_vec, ts.output.MsgQ_get_cycles_vec+ts.input.num_threads, (uint64_t)0);
            MsgQ_get_cycles_mean = MsgQ_get_cycles_sum / (double)ts.input.num_threads;
         //}
         //else if (ts.input.comp_wtime_flag == 1){
            comp_wtime_sum = accumulate(ts.output.comp_wtime_vec, ts.output.comp_wtime_vec+ts.input.num_threads, (double)0.0);
            comp_wtime_mean = comp_wtime_sum / (double)ts.input.num_threads;
         //}

         num_qGets_sum = accumulate(ts.output.num_qGets_vec, ts.output.num_qGets_vec+ts.input.num_threads, (int)0);
         num_qPuts_sum = accumulate(ts.output.num_qPuts_vec, ts.output.num_qPuts_vec+ts.input.num_threads, (int)0);
      }

      /* print output stats */
      if (verbose_output){
         printf("Solve forward-error L2-norm = %e\n"
                "Solve wall-clock time = %e\n"
                "Overall solve wtime = %e\n"
                "Setup wall-clock time = %e\n"
                "Sequential solver wall-clock time = %e\n"
                "Mean relaxations = %f\n"
                "Mean iterations = %f\n"
                "Comp wtime = %e\n"
                "MsgQ put wtime = %e\n"
                "MsgQ get wtime = %e\n"
                "MsgQ cycles = %" PRIu64 "\n"
                "Comp cycles = %" PRIu64 "\n"
                "Num qPuts = %d\n"
                "Num qGets = %d\n",
                error_x,
                ts.output.solve_wtime,
                overall_solve_wtime,
                ts.output.setup_wtime,
                seq_wtime,
                (double)num_relax_sum/(double)ts.input.num_threads,
                (double)num_iters_sum/(double)ts.input.num_threads,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
         printf("Num. level sets = %d\n"
                "Rows per level (min, max, mean) = (%d, %d, %f)\n",
                num_lvls, lvl_n_min, lvl_n_max, lvl_n_mean);
      }
      else {
         printf("%e %e %e %e %e %f %f %e %e %e %" PRIu64 " %" PRIu64 " %d %d ",
                error_x,
                ts.output.solve_wtime,
                overall_solve_wtime,
                ts.output.setup_wtime,
                seq_wtime,
                (double)num_relax_sum/(double)ts.input.num_threads,
                (double)num_iters_sum/(double)ts.input.num_threads,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
         printf("%d %d %d %f ",
                num_lvls, lvl_n_min, lvl_n_max, lvl_n_mean);
         printf("\n");
      }
   }


   char row_output_file_str[512];
   sprintf(row_output_file_str, "RowOutput_TriSolve_NumThreads%d.txt", ts.input.num_threads);
   FILE *file = fopen(row_output_file_str, "w");
   for (int i = 0; i < num_rows; i++){
      fprintf(file, "%d %e %d %d\n", 
              i, 
              ts.output.row_output[i].wtime / num_runs,
              (int)std::round((double)ts.output.row_output[i].num_spins / (double)num_runs),
              (int)std::round((double)ts.output.row_output[i].num_atomics / (double)num_runs));
   }
   fclose(file);

   free(ts.output.setup_wtime_vec);
   free(ts.output.solve_wtime_vec);
   free(ts.output.num_relax);
   free(ts.output.num_iters);

   if (ts.input.MsgQ_flag == 1){
      //if (ts.input.MsgQ_wtime_flag == 1){
         free(ts.output.MsgQ_put_wtime_vec);
         free(ts.output.MsgQ_get_wtime_vec);
      //}
      //else if (ts.input.MsgQ_cycles_flag == 1){
         free(ts.output.MsgQ_put_cycles_vec);
         free(ts.output.MsgQ_get_cycles_vec);
      //}
      //else if (ts.input.comp_wtime_flag == 1){
         free(ts.output.comp_wtime_vec);
      //}
      free(ts.output.num_qPuts_vec);
      free(ts.output.num_qGets_vec);
   }

   free(ts.output.row_output);

   return 0;
}
