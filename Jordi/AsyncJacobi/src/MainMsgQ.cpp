// General notes:
//     change double->int
//     strip down code to basics
//     num cores fixed at 16

#include <stdint.h>
//#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "mq.h"

/* CSR struct */
typedef struct{
   int *start; /* pointer to row starts */
   int *i; /* columns indices (used in COO and CSC) */
   int *j; /* row indices (used in COO and CSR) */
   int *data; /* matrix values */
   int *diag; /* diagonal elements */
   int n; /* number of rows */
   int m; /* number of columns */
   int nnz; /* number of non-zero values */
}Matrix;

int MsgQ_Put(int destination_qid, int source_data)
{
   //qPut(destination_qid, source_data)
   return 0;
}

int MsgQ_Get(int destination_qid, int *source_data)
{
   //qGet(destination_qid, &source_data);
   return 0;
}

//void PrintParMatrix(Matrix A, char *filename, int tid, int num_threads)
//{
//   for (int t = 0; t < num_threads; t++){
//      if (t == tid){
//         if (t == 0) remove(filename);
//         FILE *file = fopen(filename, "a");
//         for (int i_loc = 0; i_loc < A.n; i_loc++){
//            int i = i_loc * num_threads + tid;
//            for (int kk = A.start[i_loc]; kk < A.start[i_loc+1]; kk++){
//               int row = i+1;
//               int col = A.j[kk]+1;
//               fprintf(file, "%d %d %d\n", row, col, A.data[kk]);
//            }
//         }
//         fclose(file);
//      }
//      #pragma omp barrier
//   }
//}

void Laplace_2D_5pt(Matrix *A, /* sparse matrix data (output) */
                    int n, /* size of grid (n*n rows) */
                    int tid,
                    int num_threads)
{
   int N = n*n;
   A->n = tid < N % num_threads ? N / num_threads + 1 : N / num_threads;

   int col;

   A->m = A->n;
   int block_end = n-1;
   int block_start = 0;
   int k = 0;

   //for(int i_loc = 0; i_loc < A->n; i_loc++){
   //   int i = i_loc * num_threads + tid;
   for(int i = 0; i < N; i++){
      if (abs(i - (int)tid) % num_threads == 0){
         col = i - n;
         if (col >= 0){
            k++;
         }
         if (i > block_start){
            k++;
         }
         k++;
         if (i < block_end){
            k++;
         }
         col = i + n;
         if (col < N){
            k++;
         }
      }

      if (i == block_end){
         block_end += n;
         block_start += n;
      }
   }
   A->nnz = k;

   A->start = (int *)calloc(A->n+1, sizeof(int));
   A->j = (int *)calloc(A->nnz, sizeof(int));
   A->data = (int *)calloc(A->nnz, sizeof(int));
   A->diag = (int *)calloc(A->n, sizeof(int));

   double A_diag[A->n];

   k = 0;
   block_end = n-1;
   block_start = 0;
   A->start[0] = 0;
   //for(int i_loc = 0; i_loc < A->n; i_loc++){
   //   int i = i_loc * num_threads + tid;
   int i_loc = 0;
   for(int i = 0; i < N; i++){
      if (abs(i - tid) % num_threads == 0){
         A->diag[i_loc] = 4;
         col = i - n;
         if (col >= 0){
            A->data[k] = -1;
            A->j[k] = col;
            k++;
         }
         if (i > block_start){
            col = i - 1;
            A->data[k] = -1;
            A->j[k] = col;
            k++;
         }
         A->data[k] = 4;
         A->j[k] = i;
         k++;
         if (i < block_end){
            col = i + 1;
            A->data[k] = -1;
            A->j[k] = col;
            k++;
         }
         col = i + n;
         if (col < N){
            A->data[k] = -1;
            A->j[k] = col;
            k++;
         }
         A->start[i_loc+1] = k;
         i_loc++;
      }

      if (i == block_end){
         block_end += n;
         block_start += n;
      }
   }
}

int main (int argc, char *argv[])
{
   int num_iters = 1;
   int m = 4;

   int tid = atoi(argv[1]); // core ID
   int num_threads = 16;

   //int num_threads = atoi(argv[1]);   
   //omp_set_num_threads(num_threads);

   //#pragma omp parallel
   //{
   //   int tid = omp_get_thread_num();

      /* set up problem */
      Matrix A;
      Laplace_2D_5pt(&A, m, tid, num_threads);
      //char filename[100] = "A.txt";
      //PrintParMatrix(A, filename, tid, num_threads);
      
   
      int n_loc = A.n;
   
      int *r_loc = (int *)calloc(n_loc, sizeof(int));
      int *x_loc = (int *)calloc(n_loc, sizeof(int));
      int num_qPuts = 0, num_qGets = 0;
   
      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         int i = i_loc * num_threads + tid;
         r_loc[i_loc] = 1;
      }
   
      for (int iter = 0; iter < num_iters; iter++){
         for (int i_loc = 0; i_loc < n_loc; i_loc++){
            int i = i_loc * num_threads + tid;
            int z = r_loc[i_loc] / A.diag[i_loc];
            x_loc[i_loc] += z;
            for (int jj = A.start[i_loc]; jj < A.start[i_loc+1]; jj++){
               int j = A.j[jj];
               int y = A.data[jj] * z;
               if (i == j){
                  r_loc[i_loc] -= y;
               }
               else {
                  MsgQ_Put(j, y); // qPut()
                  num_qPuts++;
               }
            }
   
            int get_flag;
            int z_accum = 0.0, z_recv = 0.0;
            while (1){
               get_flag = MsgQ_Get(i, &z_recv); // qGet()
               if (!get_flag) break;
               z_accum += z_recv;
               num_qGets++;
            }
            r_loc[i_loc] -= z_accum;
         }
      }
   
      free(x_loc);
      free(r_loc);
   //}

   return 0;
}
