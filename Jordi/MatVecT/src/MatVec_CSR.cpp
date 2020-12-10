#include "MatVec.hpp"
#include "../../src/Matrix.hpp"

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
   #pragma omp parallel
   {
      int num_rows = A.n;
      int num_cols = A.m;
      int nnz = A.nnz;
   
      int tid = omp_get_thread_num();
      int num_threads = mv->input.num_threads;
      int i_offset = num_rows * tid;
      int j_offset = num_cols * tid;
   
      if (mv->input.AAT_flag == 1){
         if (mv->input.expand_flag == 1){
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
         if (mv->input.expand_flag == 1){
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
}
