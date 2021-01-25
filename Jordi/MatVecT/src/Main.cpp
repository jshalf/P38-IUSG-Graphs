#include "../../src/Main.hpp"
#include "MatVec.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/Misc.hpp"

int main (int argc, char *argv[])
{
   MatVecData mv;
   mv.input.num_threads = 1;
   mv.input.num_iters = 1;
   mv.input.atomic_flag = 1;
   mv.input.AAT_flag = 0;
   mv.input.expand_flag = 0;
   mv.input.coo_flag = 0;
   mv.input.MsgQ_flag = 0;
   int verbose_output = 0;
   int num_runs = 1;
   int m = 10; 
   int problem_type = PROBLEM_5PT_POISSON;
   char mat_file_str[128];

   int arg_index = 0;
   while (arg_index < argc){
      if (strcmp(argv[arg_index], "-n") == 0){
         arg_index++;
         m = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-max_iters") == 0){
         arg_index++;
         mv.input.num_iters = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-num_threads") == 0){
         arg_index++;
         mv.input.num_threads = atoi(argv[arg_index]);
      }
      else if (strcmp(argv[arg_index], "-no_atomic") == 0){
         mv.input.atomic_flag = 0;
      }
      else if (strcmp(argv[arg_index], "-num_runs") == 0){
         arg_index++;
         num_runs = std::max(1, atoi(argv[arg_index]));
      }
      else if (strcmp(argv[arg_index], "-verb_out") == 0){
         verbose_output = 1;
      }
      else if (strcmp(argv[arg_index], "-expand") == 0){
         mv.input.expand_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-AAT") == 0){
         mv.input.AAT_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-coo") == 0){
         mv.input.coo_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-MsgQ") == 0){
         mv.input.MsgQ_flag = 1;
      }
      else if (strcmp(argv[arg_index], "-problem") == 0){
         arg_index++;
         if (strcmp(argv[arg_index], "5pt") == 0){
            problem_type = PROBLEM_5PT_POISSON;
         }
         else if (strcmp(argv[arg_index], "file") == 0){
            arg_index++;
            problem_type = PROBLEM_FILE;
            strcpy(mat_file_str, argv[arg_index]);
         }
      }
      arg_index++;
   }
   
   omp_set_num_threads(mv.input.num_threads);

   CSR A;
   if (problem_type == PROBLEM_FILE){
      char A_mat_file_str[128];
      sprintf(A_mat_file_str, "%s_A.txt.bin", mat_file_str);
      freadBinaryMatrix(A_mat_file_str, &A, 1);

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
   double *y1 = (double *)calloc(num_cols, sizeof(double));
   double *y1_exact = (double *)calloc(num_cols, sizeof(double));
   double *e1 = (double *)calloc(num_cols, sizeof(double));

   double *y2, *y2_exact, *e2;
   if (mv.input.AAT_flag == 1){
      y2 = (double *)calloc(num_rows, sizeof(double));
      y2_exact = (double *)calloc(num_rows, sizeof(double));
      e2 = (double *)calloc(num_rows, sizeof(double));
   }

   if (mv.input.expand_flag == 1){
      mv.y1_expand = (double *)calloc(mv.input.num_threads * num_cols, sizeof(double));
      mv.y2_expand = (double *)calloc(mv.input.num_threads * num_rows, sizeof(double));
   }

   srand(0);
   for (int i = 0; i < num_rows; i++){
      x[i] = RandDouble(-1.0, 1.0);
   } 
   for (int run = 1; run <= num_runs; run++){
      int iter;
      double start = omp_get_wtime();
      for (int i = 0; i < num_cols; i++){
         y1[i] = 0;
         if (mv.input.AAT_flag == 1){
            y2[i] = 0;
         }
      } 
      MatVec_CSR(&mv, A, x, y1_exact);

      if (mv.input.coo_flag == 1){
         MatVecT_COO(&mv, A, x, y1, y2);
      }
      else {
         MatVecT_CSR(&mv, A, x, y1, y2);
      }

      for (int i = 0; i < num_cols; i++){
         e1[i] = y1_exact[i] - y1[i];
         //printf("%e %e\n", y1_exact[i], y1[i]);
      }
      if (mv.input.AAT_flag == 1){
         for (int i = 0; i < num_rows; i++){
            y2_exact[i] = y1_exact[i];
            e2[i] = y2_exact[i] - y2[i];
         }
      }
      mv.output.solve_wtime = omp_get_wtime() - start;
      double error1 = sqrt(InnerProd(e1, e1, num_cols))/sqrt(InnerProd(y1_exact, y1_exact, num_rows));
      double error2 = 0.0;
      if (mv.input.AAT_flag == 1) error2 = sqrt(InnerProd(e2, e2, num_rows))/sqrt(InnerProd(y2_exact, y2_exact, num_rows));
      if (verbose_output){
         printf("MatVec wall-clock time %e, AT error L2-norm %e, A error L2-norm = %e, iterations = %d\n", mv.output.solve_wtime, error1, error2, iter);
      }
      else {
         printf("%e %e %e %d\n", mv.output.solve_wtime, error1, error2, iter);
      }
   }
   return 0;

   free(x);
   free(y1);
   free(y1_exact);
   free(e1);

   if (mv.input.AAT_flag == 1){
      free(y2);
      free(y2_exact);
      free(e2);
   }

   if (mv.input.expand_flag == 1){
      free(mv.y1_expand);
      free(mv.y2_expand);
   }

   return 0;
}
