#include "MatVec.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

void MatVec_CSR(MatVecData *mv,
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

void MatVecT_CSR(MatVecData *mv,
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
   if (mv->input.MsgQ_flag == 1){
      Q.type = Q_STDQUEUE;
      q_size = 2*num_rows;
      qAlloc(&Q, q_size, NULL);
      qInitLock(&Q);
   }

   #pragma omp parallel
   {
      int tid = omp_get_thread_num();
      int num_threads = mv->input.num_threads;
      int i_offset = num_rows * tid;
      int j_offset = num_cols * tid;
   
      if (mv->input.AAT_flag == 1){
         if (mv->input.MsgQ_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  qPut(&Q, A.j[jj], z);
                  z = A.data[jj] * x[A.j[jj]];
                  qPut(&Q, num_rows+i, z);
               }
            }
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               double z;
               while(qGet(&Q, i, &z)){
                  y1[i] += z;
               }
               while(qGet(&Q, num_rows+i, &z)){
                  y2[i] += z;
               }
            }
         }
         else if (mv->input.expand_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  mv->y1_expand[j_offset + A.j[jj]] += A.data[jj] * x[i];
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
   
            #pragma omp for
            for (int i = 0; i < num_cols; i++){
               y1[i] = 0;
               for (int j = 0; j < num_threads; j++){
                  int jj = j*num_cols + i;
                  y1[i] += mv->y1_expand[jj];
                  mv->y1_expand[jj] = 0;
               }
            }
         }
         else if (mv->input.atomic_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  #pragma omp atomic
                  y1[A.j[jj]] += A.data[jj] * x[i];
                  #pragma omp atomic
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
         }
         else {
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  y1[A.j[jj]] += A.data[jj] * x[i];
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
         }
      }
      else {
         if (mv->input.MsgQ_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  qPut(&Q, A.j[jj], z);
               }
            }
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               double z;
               while(qGet(&Q, i, &z)){
                  y1[i] += z;
               }
            }
         }
         else if (mv->input.expand_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  mv->y1_expand[j_offset + A.j[jj]] += A.data[jj] * x[i];
               }
            }
   
            #pragma omp for
            for (int i = 0; i < num_cols; i++){
               y1[i] = 0;
               for (int j = 0; j < num_threads; j++){
                  int jj = j*num_cols + i;
                  y1[i] += mv->y1_expand[jj];
                  mv->y1_expand[jj] = 0;
               }
            }
         }
         else if (mv->input.atomic_flag == 1){
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  #pragma omp atomic
                  y1[A.j[jj]] += A.data[jj] * x[i];
               }
            }
         }
         else {
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.i_ptr[i]; jj < A.i_ptr[i+1]; jj++){
                  y1[A.j[jj]] += A.data[jj] * x[i];
               }
            }
         }
      }
      #pragma omp barrier
   }

   if (mv->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
