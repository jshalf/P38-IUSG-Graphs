#include "MatVec.hpp"
#include "../../src/Matrix.hpp"
#include "../../src/MsgQ.hpp"

/**************************************************************************
 * Serial y=A^Tx (combined with y=Ax if specified) where y is unknown.
 * A is in compressed sparse row (Matrix) format.
 **************************************************************************/
void MatVec_CSR(MatVecData *mv,
                Matrix A, /* sparse matrix */
                double *x, /* vector to be mulitplied with A */
                double *y /* result of Ax */
                )
{
   #pragma omp parallel
   {
      double Axi;
      int num_rows = A.n;

      #pragma omp for
      for (int i = 0; i < num_rows; i++){
         Axi = 0.0;
         for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
            Axi += A.data[jj] * x[A.j[jj]];
         }
         y[i] = Axi;
      }
   }
}

/**************************************************************************
 * Parallel y=A^Tx (combined with y=Ax if specified) where y is unknown.
 * A is in compressed sparse row (Matrix) format.
 **************************************************************************/
void MatVecT_CSR(MatVecData *mv,
                 Matrix A, /* sparse matrix */
                 double *x, /* vector to be mulitplied with A */
                 double *y1, /* result of A^Tx */
                 double *y2 /* result of Ax */
                 )
{
   int num_rows = A.n;
   int num_cols = A.m;
   int nnz = A.nnz;
   int q_size;
   Queue Q;
   int y_counts = 0;
   /* set up message queues */
   if (mv->input.MsgQ_flag == 1){
      q_size = 2*num_rows;
      qAlloc(&Q, q_size);
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
   
      if (mv->input.AAT_flag == 1){ /* compute both Ax and A^Tx  */
         if (mv->input.MsgQ_flag == 1){ /* use message queues */
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  /* compute update for element i of y1*/
                  z = A.data[jj] * x[i];
                  /* use put primitive to send update */
                  qPut(&Q, A.j[jj], z);
                  /* compute update for element i of y2*/
                  z = A.data[jj] * x[A.j[jj]];
                  /* use put primitive to send update */
                  qPut(&Q, num_rows+i, z);
               }
            }
            /* until y1 and y2 are completely updated, use get primitive to update y1 and y2 */
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
         else if (mv->input.expand_flag == 1){ /* use ``expand'' scheme */
            /* each thread accumulates a local version of y1.  y2 can be computed normally.
             * the local verions are then summed to get the desired results. */
            #pragma omp for
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  mv->y1_expand[j_offset + A.j[jj]] += A.data[jj] * x[i]; /* update local version of y1 */
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
   
            /* now sum local versions of y1 */
            #pragma omp for nowait
            for (int i = 0; i < num_cols; i++){
               y1[i] = 0;
               for (int j = 0; j < num_threads; j++){
                  int jj = j*num_cols + i;
                  y1[i] += mv->y1_expand[jj]; 
                  mv->y1_expand[jj] = 0;
               }
            }
         }
         else if (mv->input.atomic_flag == 1){ /* atomic version */
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  /* atomically accumulate element A.j[jj] of y1 */
                  #pragma omp atomic
                  y1[A.j[jj]] += A.data[jj] * x[i];
                  /* atomically accumulate element i of y2 */
                  #pragma omp atomic
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
         }
         else { /* atomic scheme with atomics removed 
                 * (for performance analysis, will not give correct result) */
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  y1[A.j[jj]] += A.data[jj] * x[i];
                  y2[i] += A.data[jj] * x[A.j[jj]];
               }
            }
         }
      }
      else { /* only compute A^Tx (for details, see comments related to computing y1 in if statement above) */
         if (mv->input.MsgQ_flag == 1){
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  double z;
                  z = A.data[jj] * x[i];
                  qPut(&Q, A.j[jj], z);
               }
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
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  mv->y1_expand[j_offset + A.j[jj]] += A.data[jj] * x[i];
               }
            }
   
            #pragma omp for nowait
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
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  #pragma omp atomic
                  y1[A.j[jj]] += A.data[jj] * x[i];
               }
            }
         }
         else {
            #pragma omp for nowait
            for (int i = 0; i < num_rows; i++){
               for (int jj = A.start[i]; jj < A.start[i+1]; jj++){
                  y1[A.j[jj]] += A.data[jj] * x[i];
               }
            }
         }
      }
   }

   if (mv->input.MsgQ_flag == 1){
      qDestroyLock(&Q);
      qFree(&Q);
   }
}
