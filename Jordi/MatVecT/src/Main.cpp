#include "../../src/Main.hpp"
#include "MatVecT.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/Misc.hpp"

int main (int argc, char *argv[])
{
   MatVecData mv;
   mv.input.num_threads = 1;
   mv.input.atomic_flag = 1;
   mv.input.AAT_flag = 0;
   mv.input.expand_flag = 0;
   mv.input.coo_flag = 0;
   mv.input.MsgQ_flag = 0;
   mv.input.comp_wtime_flag = 0;
   mv.input.MsgQ_wtime_flag = 0;
   mv.input.MsgQ_cycles_flag = 0;
   mv.input.comp_noop_flag = 0;
   mv.input.MsgQ_noop_flag = 0;
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
      else if (strcmp(argv[arg_index], "-method") == 0){ /* method */
         arg_index++;
         if (strcmp(argv[arg_index], "atomic") == 0){ /* atomic method that treats Matrix as CSC */
            mv.input.expand_flag = 0;
         }
         else if (strcmp(argv[arg_index], "expand") == 0){ /* each thread uses a local vector which are summed at the end */
            mv.input.expand_flag = 1;
         }
      }
      else if (strcmp(argv[arg_index], "-n") == 0){ /* ``size'' of matrix. n*n rows for Laplace, n rows otherwise */
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){ /* number of threads */
         arg_index++;
         mv.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){ /* use no-atomics implementations */
         mv.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){ /* number of TriSolve runs */
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){ /* verbose output */
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){ /* use message queues */
         mv.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_wtime") == 0){
         mv.input.MsgQ_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_wtime") == 0){
         mv.input.comp_wtime_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_cycles") == 0){
         mv.input.MsgQ_cycles_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ_noop") == 0){
         mv.input.MsgQ_noop_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-comp_noop") == 0){
         mv.input.comp_noop_flag = 1;
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
      printf("      expand:             baseline method that does not require atomics.\n");
      printf("-n <int value>:           size of test problem.  For 5pt, this is the length of the 2D grid, i.e., the matrix has n^2 rows.\n");
      printf("-num_threads <int value>: number of OpenMP threads.\n");
      printf("-no_atomic:               turn off atomics.  Only meant for performance measurements and will not produce a correct result.\n");
      printf("-num_runs <int value>:    number of independent runs.  Used for data collection.\n");
      printf("-verb_out:                verbose output.\n");
      printf("-MsgQ:                    use message queues instead of atomics.\n");
      printf("\n");
      return 0;
   }
   
   omp_set_num_threads(mv.input.num_threads);

   int csc_flag = 0, coo_flag = 0;
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
      Laplace_2D_5pt(mv.input, &A, m);
   }
   int num_rows = A.n;
   int num_cols = A.m;

   if (num_rows != num_cols){
      mv.input.expand_flag = 0;
   }

   double *x = (double *)calloc(num_rows, sizeof(double));
   double *y = (double *)calloc(num_cols, sizeof(double));
   double *y_exact = (double *)calloc(num_cols, sizeof(double));
   double *e = (double *)calloc(num_cols, sizeof(double));

   if (mv.input.expand_flag == 1){
      mv.y_expand = (double *)calloc(mv.input.num_threads * num_cols, sizeof(double));
   }

   mv.output.solve_wtime_vec = (double *)calloc(mv.input.num_threads, sizeof(double));
   mv.output.num_qGets_vec = (int *)calloc(mv.input.num_threads, sizeof(int));
   mv.output.num_qPuts_vec = (int *)calloc(mv.input.num_threads, sizeof(int));

   if (mv.input.MsgQ_wtime_flag == 1){
      mv.output.MsgQ_put_wtime_vec = (double *)calloc(mv.input.num_threads, sizeof(double));
      mv.output.MsgQ_get_wtime_vec = (double *)calloc(mv.input.num_threads, sizeof(double));
   }
   else if (mv.input.MsgQ_cycles_flag == 1){
      mv.output.MsgQ_put_cycles_vec = (uint64_t *)calloc(mv.input.num_threads, sizeof(uint64_t));
      mv.output.MsgQ_get_cycles_vec = (uint64_t *)calloc(mv.input.num_threads, sizeof(uint64_t));
   }
   else if (mv.input.comp_wtime_flag == 1){
      mv.output.comp_wtime_vec = (double *)calloc(mv.input.num_threads, sizeof(double));
   }

   srand(0);
   for (int i = 0; i < num_rows; i++){
      x[i] = RandDouble(-1.0, 1.0);
   } 

   MatVecT_CSR_Seq(&mv, A, x, y_exact);

   for (int run = 1; run <= num_runs; run++){
      for (int i = 0; i < num_cols; i++){
         y[i] = 0;
      }

      /* parallel MatVecT */
      MatVecT_CSR(&mv, A, x, y);

      /* compute error */
      for (int i = 0; i < num_cols; i++){
         e[i] = y_exact[i] - y[i];
         //printf("%e %e\n", y_exact[i], y[i]);
      }

      double error_norm = sqrt(InnerProd(e, e, num_cols))/sqrt(InnerProd(y_exact, y_exact, num_rows));

      double solve_wtime_sum = accumulate(mv.output.solve_wtime_vec, mv.output.solve_wtime_vec+mv.input.num_threads, 0.0);
      mv.output.solve_wtime = solve_wtime_sum / (double)mv.input.num_threads;

      double MsgQ_put_wtime_sum = 0.0, MsgQ_put_wtime_mean = 0.0;
      double MsgQ_get_wtime_sum = 0.0, MsgQ_get_wtime_mean = 0.0;
      uint64_t MsgQ_put_cycles_sum = 0;
      uint64_t MsgQ_get_cycles_sum = 0;
      double MsgQ_put_cycles_mean = 0.0;
      double MsgQ_get_cycles_mean = 0.0;
      double comp_wtime_sum = 0.0, comp_wtime_mean = 0.0;
      int num_qGets_sum = 0, num_qPuts_sum = 0;
      if (mv.input.MsgQ_flag == 1){
         if (mv.input.MsgQ_wtime_flag == 1){
            MsgQ_put_wtime_sum = accumulate(mv.output.MsgQ_put_wtime_vec, mv.output.MsgQ_put_wtime_vec+mv.input.num_threads, 0.0);
            MsgQ_put_wtime_mean = MsgQ_put_wtime_sum / (double)mv.input.num_threads;

            MsgQ_get_wtime_sum = accumulate(mv.output.MsgQ_get_wtime_vec, mv.output.MsgQ_get_wtime_vec+mv.input.num_threads, 0.0);
            MsgQ_get_wtime_mean = MsgQ_get_wtime_sum / (double)mv.input.num_threads;
         }
         else if (mv.input.MsgQ_cycles_flag == 1){
            MsgQ_put_cycles_sum = accumulate(mv.output.MsgQ_put_cycles_vec, mv.output.MsgQ_put_cycles_vec+mv.input.num_threads, 0);
            MsgQ_put_cycles_mean = MsgQ_put_cycles_sum / (double)mv.input.num_threads;
            
            MsgQ_get_cycles_sum = accumulate(mv.output.MsgQ_get_cycles_vec, mv.output.MsgQ_get_cycles_vec+mv.input.num_threads, 0);
            MsgQ_get_cycles_mean = MsgQ_get_cycles_sum / (double)mv.input.num_threads;
         }
         else if (mv.input.comp_wtime_flag == 1){
            comp_wtime_sum = accumulate(mv.output.comp_wtime_vec, mv.output.comp_wtime_vec+mv.input.num_threads, 0.0);
            comp_wtime_mean = comp_wtime_sum / (double)mv.input.num_threads;
         }
         
         num_qGets_sum = accumulate(mv.output.num_qGets_vec, mv.output.num_qGets_vec+mv.input.num_threads, 0);
         num_qPuts_sum = accumulate(mv.output.num_qPuts_vec, mv.output.num_qPuts_vec+mv.input.num_threads, 0);
      }

      /* print stats */
      if (verbose_output){
         printf("Error L2-norm %e\n"
                "MatVec wall-clock time %e\n"
                "Mean comp wtime = %e\n"
                "MsgQ put wtime = %e\n"
                "MsgQ get wtime = %e\n"
                "MsgQ put cycles = %" PRIu64 "\n"
                "MsgQ get cycles = %" PRIu64 "\n"
                "Num qPuts = %d\n"
                "Num qGets = %d\n",
                error_norm,
                mv.output.solve_wtime,
                comp_wtime_sum,
                MsgQ_put_wtime_sum,
                MsgQ_get_wtime_sum,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
      else {
         printf("%e %e %e %e %e %" PRIu64 " %" PRIu64 " %d %d\n",
                error_norm,
                mv.output.solve_wtime,
                comp_wtime_sum,
                MsgQ_put_wtime_sum,
                MsgQ_get_wtime_sum,
                MsgQ_put_cycles_sum,
                MsgQ_get_cycles_sum,
                num_qPuts_sum,
                num_qGets_sum);
      }
   }

   free(x);
   free(y);
   free(y_exact);
   free(e);

   if (mv.input.expand_flag == 1){
      free(mv.y_expand);
   }

   free(mv.output.solve_wtime_vec);
   free(mv.output.num_qGets_vec);
   free(mv.output.num_qPuts_vec);

   if (mv.input.MsgQ_wtime_flag == 1){
      free(mv.output.MsgQ_put_wtime_vec);
      free(mv.output.MsgQ_get_wtime_vec);
   }
   else if (mv.input.MsgQ_cycles_flag == 1){
      free(mv.output.MsgQ_put_cycles_vec);
      free(mv.output.MsgQ_get_cycles_vec);
   }
   else if (mv.input.comp_wtime_flag == 1){
      free(mv.output.comp_wtime_vec);
   }

   return 0;
}
