// General notes:
//     change double->int
//     strip down code to basics
//     num cores fixed at 16

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#define SOLVER_ASYNC_JACOBI 0
#define SOLVER_ASYNC_TRISOLVE 1

#define MATRIX_FULL 0
#define MATRIX_LOWER 1

#ifdef USE_OPENMP
   #include <omp.h>
   #include "../../src/MsgQ.hpp"
   #include "../../src/Misc.hpp"
   #include "../../src/Main.hpp"

/* Info about the traces,
 * used to validate Patricia's platform */
typedef struct{
   std::vector<int> phase;
   std::vector<int> source; /* source core */
   std::vector<int> dest; /* destination core */
   std::vector<size_t> msg_size; /* size of message */
   std::vector<int> data; /* message data (currently, only integers are supported) */
}TraceInfo;

/* Print the trace info */
void PrintTraces(char *filename, TraceInfo trace_loc)
{
   //remove(filename);
   FILE *file = fopen(filename, "w");
   int num_traces = trace_loc.phase.size();
   fprintf(file,    "num calls = %d\n", num_traces);
   fprintf(file, "%4s %5s %5s %5s %10s %9s", "iter", "src", "dst", "size", "data(dec)", "data(hex)\n");
   for (int i = 0; i < num_traces; i++){
      fprintf(file, "%4d %5d %5d %5zu %10d %9X\n",
                     trace_loc.phase[i],
                     trace_loc.source[i],
                     trace_loc.dest[i],
                     trace_loc.msg_size[i],
                     (int)trace_loc.data[i],
                     (int)trace_loc.data[i]);
   }
   fclose(file);
}
#else
   #include "mq.h"
#endif

int num_threads = 16;

/* CSR struct */
typedef struct{
   int *start; /* pointer to row starts - Size=16 integers*/
   int *i;     /* columns indices (used in COO and CSC)*/
   int *j;     /* row indices (used in COO and CSR) - Size= 64 integers*/
   int *data;  /* matrix values - Size=64 integers*/
   int *diag;  /* diagonal elements - Size=16 integers*/
   int n;      /* number of rows */
   int m;      /* number of columns */
   int nnz;    /* number of non-zero values */
}Matrix;

/* Message struct */
typedef struct {
   int src_id; /* Source core */
   uint32_t data; /* Message data */
}Message;

#ifdef USE_OPENMP
MessageQueue<Message> *Q;
#endif

/* qPut wrapper that calls software or hardware message queues */
int MsgQ_Put(int destination_qid, uint32_t source_data)
{
#ifdef USE_OPENMP
   int tid = omp_get_thread_num();
   Message msg{tid, source_data};
   Q->qPut(destination_qid, msg);
   return 0;
#else
   /*qPut((uint32_t)destination_qid, (uint32_t)source_data);*/
   qPut(destination_qid, source_data);
   return 0;
#endif
}

/* qGet wrapper that calls software or hardware message queues */
int MsgQ_Get(int destination_qid, uint32_t *source_data)
{
#ifdef USE_OPENMP
   Message msg;
   int received_data_poll;
   /* Check if queue is empty */
   received_data_poll = Q->qGet(destination_qid, &msg);
   if (received_data_poll == 1){
      /* If queue has data, copy popped data */
      *source_data = msg.data;
      return msg.src_id;
   }
   return -1;
#else
   uint32_t temp;
   uint32_t received_data_poll;
   uint32_t source_id;
   uint32_t source_mask;
   source_mask = 0X00FC0000;
   qPoll((uint32_t)destination_qid, received_data_poll);
   /* Check if queue is empty */
   if (received_data_poll == 1){
      return -1;
   }
   else {
      /* If queue has data, get source core using a mask */
      source_id = received_data_poll & source_mask;
      source_id = source_id >> 17;
      /* Pop from queue */
      qGet((uint32_t)destination_qid, temp);
      /* Copy poppwed data */
      *source_data = temp;
      return 1;
   }
#endif
}

/* Function that maps sparse matrix row indices to cores */
int row_to_core(int i)
{
   /* Round-robin distribution of rows to cores */
   return i % num_threads;
}

/* Function that maps local to global sparse matrix row indices */
int loc_to_glob_row(int i_loc, int tid)
{
   /* Round-robin distribution of rows to cores */
   return i_loc * num_threads + tid;
}

#ifdef USE_OPENMP
/* Print sparse matrix (used for debugging) */
void PrintParMatrix(Matrix A, char *filename, int tid, int num_threads)
{
   for (int t = 0; t < num_threads; t++){
      if (t == tid){
         if (t == 0) remove(filename);
         FILE *file = fopen(filename, "a");
         for (int i_loc = 0; i_loc < A.n; i_loc++){
            int i = i_loc * num_threads + tid;
            for (int kk = A.start[i_loc]; kk < A.start[i_loc+1]; kk++){
               int row = i+1;
               int col = A.j[kk]+1;
               fprintf(file, "%d %d %d\n", row, col, A.data[kk]);
            }
         }
         fclose(file);
      }
      #pragma omp barrier
   }
}
#endif

int main (int argc, char *argv[])
{
   int num_iters = 1;
   int m = 4;
   int n = m * m;
   int async_flag = 1;
   int solver_type = SOLVER_ASYNC_TRISOLVE;
   int matrix_type = MATRIX_FULL;

#ifdef USE_OPENMP
   vector<int> x(n);
   /* Set number of OpenMP threads */
   omp_set_num_threads(num_threads);
   /* Initialize software message queues */
   Q = new MessageQueue<Message>(num_threads);
   #pragma omp parallel
   {
      /* Thread id */
      int tid = omp_get_thread_num();
#else
      /* Core id */
      int tid = atoi(argv[1]);
#endif

/******************************************************
 *                     Setup matrix                   *
 ******************************************************/
      Matrix A;
      //Laplace_2D_5pt(&A, m, tid, num_threads);
      //char filename[100] = "A.txt";
      //PrintParMatrix(A, filename, tid, num_threads);

      
      int N = m*m;
      A.n = tid < N % num_threads ? N / num_threads + 1 : N / num_threads;
   
      int col;
   
      A.m = A.n;
      int block_end = m-1;
      int block_start = 0;
      int k = 0;
  
      /* Construct five-point discretization of the 2D Laplace equation */ 
      //for(int i_loc = 0; i_loc < A.n; i_loc++){
      //   int i = i_loc * num_threads + tid;
      for(int i = 0; i < N; i++){
         if (row_to_core(i) == tid){
            col = i - m;
            if (col >= 0){
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }
      
               if (insert_nnz){
                  k++;
               }
            }
            if (i > block_start){
               col = i - 1;
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  k++;
               }
            }
            k++;
            if (i < block_end){
               col = i + 1;
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }
               
               if (insert_nnz){
                  k++;
               }
            }
            col = i + m;
            if (col < N){
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  k++;
               }
            }
         }
   
         if (i == block_end){
            block_end += m;
            block_start += m;
         }
      }
      A.nnz = k;
   
      int A_start[A.n+1];
      int A_j[A.nnz];
      int A_data[A.nnz];
   
      k = 0;
      block_end = m-1;
      block_start = 0;
      A_start[0] = 0;

      //for(int i_loc = 0; i_loc < A.n; i_loc++){
      //   int i = i_loc * num_threads + tid;
      int i_loc = 0;
      for(int i = 0; i < N; i++){
         if (row_to_core(i) == tid){
            //A.diag[i_loc] = 4;
            col = i - m;
            if (col >= 0){
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  A_data[k] = -1;
                  A_j[k] = col;
                  k++;
               }
            }
            if (i > block_start){
               col = i - 1;
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  A_data[k] = -1;
                  A_j[k] = col;
                  k++;
               }
            }
            A_data[k] = 1;
            A_j[k] = i;
            k++;
            if (i < block_end){
               col = i + 1;
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  A_data[k] = -1;
                  A_j[k] = col;
                  k++;
               }
            }
            col = i + m;
            if (col < N){
               bool insert_nnz = false;
               if (matrix_type == MATRIX_LOWER){
                  if (col < i) insert_nnz = true;
               }
               else {
                  insert_nnz = true;
               }

               if (insert_nnz){
                  A_data[k] = -1;
                  A_j[k] = col;
                  k++;
               }
            }
            A_start[i_loc+1] = k;
            i_loc++;
         }
   
         if (i == block_end){
            block_end += m;
            block_start += m;
         }
      }

/***********************
 * end of matrix setup 
 ***********************/
      
      int n_loc = A.n;
      int r_loc[n_loc];
      int x_loc[n_loc];
      int num_qPuts = 0, num_qGets = 0;
      int num_recvs = 0;

      /* Set initial approximation x (Jacobi only), right-had side, and residual (JJacobi only) */
      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         for (int jj = A_start[i_loc]; jj < A_start[i_loc+1]; jj++){
            int j = A_j[jj];
            if (row_to_core(j) != tid){
               num_recvs++;
            }
         }
         int i = loc_to_glob_row(i_loc, tid);
         if (solver_type == SOLVER_ASYNC_TRISOLVE){
            r_loc[i_loc] = i;
            x_loc[i_loc] = r_loc[i_loc];
         }
         else {
            r_loc[i_loc] = i;
            x_loc[i_loc] = 0;
         }
      }


#ifdef USE_OPENMP
      int print_traces = 1;
      TraceInfo trace_put, trace_get;
#endif

/************
 * Trisolve
 ************/
      if (solver_type == SOLVER_ASYNC_TRISOLVE){
         for (int i_loc = 0; i_loc < n_loc; i_loc++){ /* loop over rows */
            int i = loc_to_glob_row(i_loc, tid);
            int jj_start = A_start[i_loc];
            int jj_end = A_start[i_loc+1];
            int jj = jj_start;
            int row_count_i = 0;
            /* Count non-zeros in lower part */
            for (jj = jj_start; jj < jj_end; jj++){
               if (A_j[jj] >= i) break;
               row_count_i++;
            }
            jj = jj_start;
            /* Loop until all dependencies have been satisfied */
            while (row_count_i > 0){
               for (jj = jj_start; jj < jj_end; jj++){
                  /* Stop if we are in upper part of matrix */
                  if (A_j[jj] >= i) break;
                  /* If message has arrived, update x_loc[i_loc] */
                  uint32_t xj;
                  int src_tid = MsgQ_Get(tid, &xj); // qGet()
                  if (src_tid > -1) {
                     x_loc[i_loc] -= A_data[jj] * (int)xj;
                     row_count_i--;
                     num_qGets++;
#ifdef USE_OPENMP
                     if (print_traces == 1){
                        trace_get.phase.push_back(0);
                        trace_get.source.push_back(src_tid);
                        trace_get.dest.push_back(tid);
                        trace_get.msg_size.push_back(sizeof(int));
                        trace_get.data.push_back((int)xj);
                     }
#endif
                  }
               }
            }
            x_loc[i_loc] /= A_data[jj];
            jj++;

            /* Send x_loc[i_loc] to dependencies */
            int xi = x_loc[i_loc];
            for (; jj < jj_end; jj++){
               int dest_tid = row_to_core(A_j[jj]);
               MsgQ_Put(dest_tid, xi);
               num_qPuts++;
#ifdef USE_OPENMP
               if (print_traces == 1){
                  trace_put.phase.push_back(0);
                  trace_put.source.push_back(tid);
                  trace_put.dest.push_back(dest_tid);
                  trace_put.msg_size.push_back(sizeof(int));
                  trace_put.data.push_back(xi);
               }
#endif
            }
         }
      }
/****************************************************************
 * Jacobi (currently doesn't converge, we need double precision)
 ****************************************************************/    
      else {
         for (int iter = 0; iter < num_iters; iter++){
            for (int i_loc = 0; i_loc < n_loc; i_loc++){
               int i = loc_to_glob_row(i_loc, tid);
               int z = r_loc[i_loc];// / A.diag[i_loc];
               x_loc[i_loc] += z;
               for (int jj = A_start[i_loc]; jj < A_start[i_loc+1]; jj++){
                  int j = A_j[jj];
                  int y = A_data[jj] * z;
                  if (i == j){
                     r_loc[i_loc] -= y;
                  }
                  else {
                     MsgQ_Put(j, y);
                     num_qPuts++;
#ifdef USE_OPENMP
                     if (print_traces == 1){
                        trace_put.phase.push_back(iter);
                        trace_put.source.push_back(tid);
                        trace_put.dest.push_back(j);
                        trace_put.msg_size.push_back(sizeof(int));
                        trace_put.data.push_back((int)y);
                     }
#endif
                  }
               }
            }
   
            int z_accum = 0.0;
            uint32_t z_recv = 0.0;
            if (async_flag == 1){
               while (1){
                  int src_tid = MsgQ_Get(tid, &z_recv);
                  if (src_tid == -1) break;
                  z_accum += (int)z_recv;
                  num_qGets++;
#ifdef USE_OPENMP
                  if (print_traces == 1){
                     trace_get.phase.push_back(iter);
                     trace_get.source.push_back(src_tid);
                     trace_get.dest.push_back(tid);
                     trace_get.msg_size.push_back(sizeof(int));
                     trace_get.data.push_back((int)z_recv);
                  }
#endif
               }
               r_loc[0] -= z_accum;
            }
            /* Synchronous Jacobi only works for one iteration,
               we need to implement a barrier for correctness */ 
            else {
               int recv_count = num_recvs;
   
               while (recv_count){
                  int src_tid = MsgQ_Get(tid, &z_recv); 
                  if (src_tid > -1) {
                     recv_count--;
                     z_accum += (int)z_recv;
                     num_qGets++;
#ifdef USE_OPENMP
                     if (print_traces == 1){
                        trace_get.phase.push_back(iter);
                        trace_get.source.push_back(src_tid);
                        trace_get.dest.push_back(tid);
                        trace_get.msg_size.push_back(sizeof(int));
                        trace_get.data.push_back((int)z_recv);
                     }
#endif
                  }
               }
               r_loc[0] -= z_accum;
            }
         }
#ifdef USE_OPENMP
      }
      for (int i_loc = 0; i_loc < n_loc; i_loc++){
         int i = loc_to_glob_row(i_loc, tid);
         x[i] = x_loc[i_loc];
      }

      if (print_traces == 1){
         char trace_filename[100];
         sprintf(trace_filename, "trace_files/put_traces_%d", tid);
         PrintTraces(trace_filename, trace_put);
         sprintf(trace_filename, "trace_files/get_traces_%d", tid);
         PrintTraces(trace_filename, trace_get);
      }
   }
   delete Q;
#else
      }
#endif

   return 0;
}
