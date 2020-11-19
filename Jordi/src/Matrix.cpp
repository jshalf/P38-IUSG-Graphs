#include "Matrix.hpp"

void Laplace_2D_5pt(InputData input, CSR *A, int n)
{
   int col;
   int N = n*n;

   A->n = N;
   A->m = N;
   A->nnz = 5*n*n - 4*n;
   A->i_ptr = (int *)calloc(N+1, sizeof(int));
   A->j = (int *)calloc(A->nnz, sizeof(int));
   A->data = (double *)calloc(A->nnz, sizeof(double));
   A->diag = (double *)calloc(A->n, sizeof(double));

   int block_end = n-1;
   int block_start = 0;
   int k = 0;
   A->i_ptr[0] = 0;
   for(int i = 0; i < N; i++){
      A->diag[i] = 4.0;
      col = i - n;
      if (col >= 0){
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      if (i > block_start){
         col = i - 1;
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      A->data[k] = 4.0;
      A->j[k] = i;
      k++;
      if (i < block_end){
         col = i + 1;
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      col = i + n;
      if (col < N){
         A->data[k] = -1.0;
         A->j[k] = col;
         k++;
      }
      A->i_ptr[i+1] = k;

      if (i == block_end){
         block_end += n;
         block_start += n;
      }
   }

   if (input.coo_flag == 1){
      A->i = (int *)calloc(A->nnz, sizeof(int));
      k = 0;
      for (int i = 0; i < A->n; i++){
         for (int jj = A->i_ptr[i]; jj < A->i_ptr[i+1]; jj++){
            A->i[k] = i;
            k++;
         }
      }
   }
}
