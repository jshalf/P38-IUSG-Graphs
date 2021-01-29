#include "MatVec.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

void MatVec_COO(MatVecData *mv,
                CSR A,
                double *x,
                double *y)
{
   #pragma omp parallel
   {
      double Axi;
      int num_rows = A.n;

      #pragma omp for
      for (int i = 0; i < num_rows; i++){
         Axi = 0.0;
         for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
            Axi += A.data[jj] * x[A.j[jj]];
         }
         y[i] = Axi;
      }
   }
}

void MatVecT_COO(MatVecData *mv,
                 CSR A,
                 double *x,
                 double *y1,
                 double *y2)
{
   int num_rows = A.n;
   int num_cols = A.m;
   int nnz = A.nnz;
   int q_size;
   Queue Q;
   int y_counts = 0;
   if (mv->input.MsgQ_flag == 1){
      Q.type = Q_STDQUEUE;
      q_size = 2*num_rows;
      qAlloc(&Q, q_size, NULL);
      qInitLock(&Q);
      if (mv->input.AAT_flag == 1){
         y_counts = 2 * A.nnz;
      }
      else {
         y_counts = A.nnz;
      }
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = mv->input.num_threads;
      int i_offset = num_rows * tid;
      int j_offset = num_cols * tid;
   
      if (mv->input.AAT_flag == 1){
         if (mv->input.MsgQ_flag == 1){
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               double z;
               z = A.data[k] * x[A.i[k]];
               qPut(&Q, A.j[k], z);
               z = A.data[k] * x[A.j[k]];
               qPut(&Q, num_rows+A.i[k], z);
            }
            while (y_counts > 0){
               #pragma omp for nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  while(qGet(&Q, i, &z)){
                     y1[i] += z;
                     y_counts--;
                  }
                  while(qGet(&Q, num_rows+i, &z)){
                     y2[i] += z;
                     y_counts--;
                  }
               }
            }
            #pragma omp barrier
         }
         else if (mv->input.expand_flag == 1){
            #pragma omp for
            for (int k = 0; k < nnz; k++){
               mv->y1_expand[j_offset + A.j[k]] += A.data[k] * x[A.i[k]];
               mv->y2_expand[i_offset + A.i[k]] += A.data[k] * x[A.j[k]];
            }
   
            #pragma omp for nowait
            for (int k = 0; k < num_cols; k++){
               y1[k] = 0;
               for (int kk = 0; kk < num_threads; kk++){
                  int jj = kk*num_cols + k;
                  y1[k] += mv->y1_expand[jj];
                  mv->y1_expand[jj] = 0;
               }
            }
   
            #pragma omp for nowait
            for (int k = 0; k < num_rows; k++){
               y2[k] = 0;
               for (int kk = 0; kk < num_threads; kk++){
                  int jj = kk*num_rows + k;
                  y2[k] += mv->y2_expand[jj];
                  mv->y2_expand[jj] = 0;
               }
            }
         }
         else if (mv->input.atomic_flag == 1){
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               #pragma omp atomic
               y1[A.j[k]] += A.data[k] * x[A.i[k]];
               #pragma omp atomic
               y2[A.i[k]] += A.data[k] * x[A.j[k]];
            }
         }
         else {
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               y1[A.j[k]] += A.data[k] * x[A.i[k]];
               y2[A.i[k]] += A.data[k] * x[A.j[k]];
            }
         }
      }
      else {
         if (mv->input.MsgQ_flag == 1){
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               double z;
               z = A.data[k] * x[A.i[k]];
               qPut(&Q, A.j[k], z);
            }
            while (y_counts > 0){
               #pragma omp for nowait
               for (int i = 0; i < num_rows; i++){
                  double z;
                  while(qGet(&Q, i, &z)){
                     y1[i] += z;
                     y_counts--;
                  }
               }
            }
            #pragma omp barrier
         }
         else if (mv->input.expand_flag == 1){
            #pragma omp for
            for (int k = 0; k < nnz; k++){
               mv->y1_expand[j_offset + A.j[k]] += A.data[k] * x[A.i[k]];
            }
            #pragma omp for nowait
            for (int k = 0; k < num_cols; k++){
               y1[k] = 0;
               for (int kk = 0; kk < num_threads; kk++){
                  int jj = kk*num_cols + k;
                  y1[k] += mv->y1_expand[jj];
                  mv->y1_expand[jj] = 0;
               }
            }
         }
         else if (mv->input.atomic_flag == 1){
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               #pragma omp atomic
               y1[A.j[k]] += A.data[k] * x[A.i[k]];
            }
         }
         else {
            #pragma omp for nowait
            for (int k = 0; k < nnz; k++){
               y1[A.j[k]] += A.data[k] * x[A.i[k]];
            }
         }
      }
   }

   if (mv->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
