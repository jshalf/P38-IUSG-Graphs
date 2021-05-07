#include "Main.hpp"
#include "MsgQ.hpp"
#include "Misc.hpp"

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
   omp_set_lock(&(Q->lock[destinationQID])); /* acquire lock */

   Q->q[destinationQID].push(sourceData); /* push to queue */

   omp_unset_lock(&(Q->lock[destinationQID])); /* release lock */
}

int qPoll(Queue *Q,
          int destinationQID, /* destination queue id */
          double *destinationData /* destination data */
          )
{
   int flag = 0;
   double x;
   //omp_set_lock(&(Q->lock[i])); /* acquire lock */

   if (omp_test_lock(&(Q->lock[destinationQID]))){
      if (!(Q->q[destinationQID].empty())){ /* if the queue is empty, return zero flag */
         /* if queue is not empty, get front and flag of one */
         *destinationData = Q->q[destinationQID].front();
         flag = 1;
      }

      omp_unset_lock(&(Q->lock[destinationQID])); /* release lock */
   }

   return flag;
}

void qWait(Queue *Q,
           int destinationQID, /* destination queue id */
           double *destinationData /* destination data */
           )
{
   double x;

   while (qPoll(Q, destinationQID, &x) == 0); /* poll until queue is not empty */

   *destinationData = x;
}

/* get primitive */
int qGet(Queue *Q,
         int destinationQID, /* destination queue id */
         double *destinationData /* destination data */
         )
{
   double x;
   int flag = 0;
   //omp_set_lock(&(Q->lock[i])); /* acquire lock */

   if (omp_test_lock(&(Q->lock[destinationQID]))){
      if (!(Q->q[destinationQID].empty())){ /* if the queue is empty, return zero flag */
         /* if queue is not empty, pop front and flag of one */
         *destinationData = Q->q[destinationQID].front();
         Q->q[destinationQID].pop(); /* pop front */
         flag = 1;
      }

      omp_unset_lock(&(Q->lock[destinationQID])); /* release lock */
   }

   return flag;
}
