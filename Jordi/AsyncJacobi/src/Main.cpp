#include "../../src/Main.hpp"
#include "../../src/Misc.hpp"
#include "../../src/Matrix.hpp"
#include "Jacobi.hpp"

int main (int argc, char *argv[])
{
   /* set defaults */
   SolverData solver;
   solver.input.solver_type = ASYNC_JACOBI;
   solver.input.num_threads = 1;
   solver.input.num_iters = 200;
   solver.input.atomic_flag = 1;
   solver.input.MsgQ_flag = 0;
   solver.input.mat_storage_type = MATRIX_STORAGE_CSR;
   solver.input.comp_wtime_flag = 0;
   solver.input.MsgQ_wtime_flag = 0;
   solver.input.MsgQ_cycles_flag = 0;
   solver.input.comp_noop_flag = 0;
   solver.input.MsgQ_noop_flag = 0;
   solver.input.symm_flag = 1;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   double w = 1.0;
   int problem_type = PROBLEM_5PT_POISSON;
   char mat_file_str[128];

   int arg_index = 0;
   int print_usage = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-problem") == 0){ /* test problem */
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){ /* five-point centered-difference Laplace problem */
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "file") == 0){ /* read matrix from binary file */
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      else if (strcmp(argv[arg_index], "-solver") == 0){ /* solver name */
         arg_index++;
         if (strcmp(argv[arg_index], "sj") == 0){ /* synchronous Jacobi */
            solver.input.solver_type = SYNC_JACOBI;
         }
         else if (strcmp(argv[arg_index], "aj") == 0){ /* asynchronous Jacobi */
            solver.input.solver_type = ASYNC_JACOBI;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* size of matrix. n*n rows for Laplace, n rows otherwise */
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_iters") == 0){ /* max number of iterations */
         arg_index++;
         solver.input.num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         solver.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-w") == 0){ /* weight for weighsolver Jacobi (currently not used) */
         arg_index++;
         w = atof(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use atomics for async Jacobi */
         solver.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of separate Jacobi runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues in async solvers */
         solver.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-help") == 0){ /* print command line options */
         print_usage = 1;
      }
      else if (strcmp(argv[arg_index], "-sp_store_type") == 0){ /* matrix storage type */
         arg_index++;
         if (strcmp(argv[arg_index], "csr") == 0){ /* compressed sparse row */
            solver.input.mat_storage_type = MATRIX_STORAGE_CSR;
         }
         else if (strcmp(argv[arg_index], "csc") == 0){ /* compressed sparse column */
            solver.input.mat_storage_type = MATRIX_STORAGE_CSC;
         }
      }
      else if (strcmp(argv[arg_index], "-MsgQ_wtime") == 0){
         solver.input.MsgQ_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_wtime") == 0){
         solver.input.comp_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_cycles") == 0){
         solver.input.MsgQ_cycles_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_noop") == 0){
         solver.input.MsgQ_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_noop") == 0){
         solver.input.comp_noop_flag = 1;
      }
      arg_index++;
   }

   if (print_usage == 1){
      printf("\n");
      printf("-problem <problem_name>:  test problem.\n");
      printf("      5pt:                five-point centered difference discretization of the Poisson equation.\n");
      printf("      file:               read matrix from binary file.\n"
             "                          Must be in (i,j,val) format starting with (num rows, num cols, num nnz) on first line.\n");
      printf("-solver <solver_name>:    method for solving Ax=b.\n");
      printf("      sj:                 synchronous Jacobi.\n");
      printf("      aj:                 asynchronous Jacobi.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measuremensolver and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("-num_iters <int value>:   number of Jacobi iterations.\n");
      printf("\n");
      return 0;
   }

   if (solver.input.solver_type == SYNC_JACOBI){
      solver.input.MsgQ_flag = 0;
   }
   
   omp_set_num_threads(solver.input.num_threads);

   int csc_flag = 0, coo_flag = 0;
   if (solver.input.mat_storage_type == MATRIX_STORAGE_CSC){
      csc_flag = 1;
   }
   int include_diag = 1;

   /* set up problem */
   Matrix A;
   if (problem_type == PROBLEM_FILE){
      char A_mat_file_str[128];
      sprintf(A_mat_file_str, "%s", mat_file_str);
      freadBinaryMatrix(A_mat_file_str, &A, include_diag, csc_flag, coo_flag, MATRIX_NONSYMMETRIC);
      //char A_outfile[128];
      //sprintf(A_outfile, "./matlab/A.txt");
      //PrintCOO(A, A_outfile, 0);
   }
   else {
      Laplace_2D_5pt(solver.input, &A, m);
   }

   int n = A.n;
   double *x = (double *)calloc(n, sizeof(double));
   double *b = (double *)calloc(n, sizeof(double));

   solver.output.solve_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
   solver.output.num_qGets_vec = (int *)calloc(solver.input.num_threads, sizeof(int));
   solver.output.num_qPuts_vec = (int *)calloc(solver.input.num_threads, sizeof(int));

   if (solver.input.MsgQ_flag == 1){
      //if (solver.input.MsgQ_wtime_flag == 1){
         solver.output.MsgQ_put_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
         solver.output.MsgQ_get_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
      //}
      //else if (solver.input.MsgQ_cycles_flag == 1){
         solver.output.MsgQ_put_cycles_vec = (uint64_t *)calloc(solver.input.num_threads, sizeof(uint64_t));
         solver.output.MsgQ_get_cycles_vec = (uint64_t *)calloc(solver.input.num_threads, sizeof(uint64_t));
      //}
      //else if (solver.input.comp_wtime_flag == 1){
         solver.output.comp_wtime_vec = (double *)calloc(solver.input.num_threads, sizeof(double));
      //}
   }

   for (int run = 1; run <= num_runs; run++){
      for (int t = 0; t < solver.input.num_threads; t++){
         solver.output.solve_wtime_vec[t] = 0.0;
      }
      srand(0);
      for (int i = 0; i < n; i++){
         b[i] = RandDouble(-1.0, 1.0);
         x[i] = 0;
      }
      /* run Jacobi solver */
      Jacobi(&solver, A, b, &x);

      double solve_wtime_sum = accumulate(solver.output.solve_wtime_vec, solver.output.solve_wtime_vec+solver.input.num_threads, (double)0.0);
      solver.output.solve_wtime = solve_wtime_sum / (double)solver.input.num_threads;

      double res_norm;
      if (solver.input.mat_storage_type == MATRIX_STORAGE_CSC){
         res_norm = Residual2Norm_CSC(A, x, b);
      }
      else {
         res_norm = Residual2Norm(A, x, b);
      }

      double MsgQ_put_wtime_sum = 0.0, MsgQ_put_wtime_mean = 0.0;
      double MsgQ_get_wtime_sum = 0.0, MsgQ_get_wtime_mean = 0.0;
      uint64_t uint64_t_accum_start = 0;
      uint64_t MsgQ_put_cycles_sum = 0;
      uint64_t MsgQ_get_cycles_sum = 0;
      double MsgQ_put_cycles_mean = 0.0;
      double MsgQ_get_cycles_mean = 0.0;
      double comp_wtime_sum = 0.0, comp_wtime_mean = 0.0;
      int num_qGets_sum = 0, num_qPuts_sum = 0;
      if (solver.input.MsgQ_flag == 1){
         //if (solver.input.MsgQ_wtime_flag == 1){
            MsgQ_put_wtime_sum = accumulate(solver.output.MsgQ_put_wtime_vec, solver.output.MsgQ_put_wtime_vec+solver.input.num_threads, (double)0.0);
            MsgQ_put_wtime_mean = MsgQ_put_wtime_sum / (double)solver.input.num_threads;

            MsgQ_get_wtime_sum = accumulate(solver.output.MsgQ_get_wtime_vec, solver.output.MsgQ_get_wtime_vec+solver.input.num_threads, (double)0.0);
            MsgQ_get_wtime_mean = MsgQ_get_wtime_sum / (double)solver.input.num_threads;
         //}
         //else if (solver.input.MsgQ_cycles_flag == 1){
            MsgQ_put_cycles_sum = accumulate(solver.output.MsgQ_put_cycles_vec, solver.output.MsgQ_put_cycles_vec+solver.input.num_threads, (uint64_t)0);
            MsgQ_put_cycles_mean = MsgQ_put_cycles_sum / (double)solver.input.num_threads;

            MsgQ_get_cycles_sum = accumulate(solver.output.MsgQ_get_cycles_vec, solver.output.MsgQ_get_cycles_vec+solver.input.num_threads, (uint64_t)0);
            MsgQ_get_cycles_mean = MsgQ_get_cycles_sum / (double)solver.input.num_threads;
         //}
         //else if (solver.input.comp_wtime_flag == 1){
            comp_wtime_sum = accumulate(solver.output.comp_wtime_vec, solver.output.comp_wtime_vec+solver.input.num_threads, (double)0.0);
            comp_wtime_mean = comp_wtime_sum / (double)solver.input.num_threads;
         //}

         num_qGets_sum = accumulate(solver.output.num_qGets_vec, solver.output.num_qGets_vec+solver.input.num_threads, (int)0);
         num_qPuts_sum = accumulate(solver.output.num_qPuts_vec, solver.output.num_qPuts_vec+solver.input.num_threads, (int)0);
      }

      /* print solver stats */
      if (verbose_output){
         printf("Rel res. 2-norm = %e\n"
                "Solve wall-clock time = %e\n"
                "Mean comp wtime = %e\n"
                "MsgQ put wtime = %e\n"
                "MsgQ get wtime = %e\n"
                "MsgQ put cycles = %" PRIu64 "\n"
                "MsgQ get cycles = %" PRIu64 "\n"
                "Num qPuts = %d\n"
                "Num qGets = %d\n",
                res_norm,
                solver.output.solve_wtime,
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
                res_norm,
                solver.output.solve_wtime,
                comp_wtime_mean,
                MsgQ_put_wtime_mean,
                MsgQ_get_wtime_mean,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
   }

   free(x);
   free(b);

   free(solver.output.solve_wtime_vec);
   free(solver.output.num_qGets_vec);
   free(solver.output.num_qPuts_vec);

   if (solver.input.MsgQ_flag == 1){
      //if (solver.input.MsgQ_wtime_flag == 1){
         free(solver.output.MsgQ_put_wtime_vec);
         free(solver.output.MsgQ_get_wtime_vec);
      //}
      //else if (solver.input.MsgQ_cycles_flag == 1){
         free(solver.output.MsgQ_put_cycles_vec);
         free(solver.output.MsgQ_get_cycles_vec);
      //}
      //else if (solver.input.comp_wtime_flag == 1){
         free(solver.output.comp_wtime_vec);
      //}
   }

   return 0;
}
