#include "Main.hpp"
#include "MsgQ.hpp"

/* initialize mutex lock */
void qInitLock(Queue *Q)
{
   for (int i = 0; i < Q->size; i++){
      omp_init_lock(&(Q->lock[i]));
   }
}

/* destroy mutex lock */
void qDestroyLock(Queue *Q)
{
   for (int i = 0; i < Q->size; i++){
      omp_destroy_lock(&(Q->lock[i]));
   }
}

/* initialize queue data */
void qAlloc(Queue *Q,
            int n /* length of vector of queues */
            )
{
   Q->size = n;
   Q->q = new queue<double>[Q->size];
   Q->lock = (omp_lock_t *)malloc(Q->size * sizeof(omp_lock_t));
}

/* free queue data */
void qFree(Queue *Q)
{
   delete [] Q->q;
}

/* put primitive */
void qPut(Queue *Q, 
          int destinationQID, /* destination queue id */
          double sourceData /* source data */
          )
{
   int i = destinationQID;

   omp_set_lock(&(Q->lock[i])); /* acquire lock */

   Q->q[i].push(sourceData); /* push to queue */

   omp_unset_lock(&(Q->lock[i])); /* release lock */
}

int qPoll(Queue *Q,
          int destinationQID, /* destination queue id */
          double *destinationData /* destination data */
          )
{
   int flag = 0;
   int i = destinationQID;
   double x;
    
   //omp_set_lock(&(Q->lock[i])); /* acquire lock */

   if (omp_test_lock(&(Q->lock[i]))){
      if (!(Q->q[i].empty())){ /* if the queue is empty, return zero flag */
         /* if queue is not empty, get front and flag of one */
         *destinationData = Q->q[i].front();
         flag = 1;
      }

      omp_unset_lock(&(Q->lock[i])); /* release lock */
   }

   return flag;
}

void qWait(Queue *Q,
           int destinationQID, /* destination queue id */
           double *destinationData /* destination data */
           )
{
   int q_empty = 1;
   int i = destinationQID;
   double x;

   while (qPoll(Q, i, &x) == 0); /* poll until queue is not empty */

   *destinationData = x;
}

/* get primitive */
int qGet(Queue *Q,
         int destinationQID, /* destination queue id */
         double *destinationData /* destination data */
         )
{
   int q_empty = 1;
   int i = destinationQID;
   double x;
   int flag = 0;

   //omp_set_lock(&(Q->lock[i])); /* acquire lock */

   if (omp_test_lock(&(Q->lock[i]))){
      if (!(Q->q[i].empty())){ /* if the queue is empty, return zero flag */
         /* if queue is not empty, pop front and flag of one */
         *destinationData = Q->q[i].front();
         Q->q[i].pop(); /* pop front */
         flag = 1;
      }

      omp_unset_lock(&(Q->lock[i])); /* release lock */
   }

   return flag;
}

/* accumulation (assumes single-write/mulitple-read) */
void qAccum(Queue *Q,
            int destinationQID, /* destination queue id */
            double y /* accumulation data */
            )
{
   double x;
   qWait(Q, destinationQID, &x); /* get the queue data */
   x += y; /* accumulate */
   qPut(Q, destinationQID, x); /* push the accum result to the queue */
   qGet(Q, destinationQID, &x); /* pop off old data */
}
