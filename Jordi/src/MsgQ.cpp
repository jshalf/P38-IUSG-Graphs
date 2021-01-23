#include "Main.hpp"
#include "MsgQ.hpp"

void qPopArray(Queue *Q, int destinationQID)
{
   int i = destinationQID;

   for (int j = 0; j < Q->position[i]; j++){
      Q->a[i][j] = Q->a[i][j+1];
   }
}

void qInitLock(Queue *Q)
{
   omp_init_lock(&(Q->lock));
}

void qDestroyLock(Queue *Q)
{
   omp_destroy_lock(&(Q->lock));
}

void qAlloc(Queue *Q, int n, int *q_lens)
{
   Q->size = n;
   if (Q->type == Q_ARRAY){ 
      Q->position = (int *)malloc(Q->size * sizeof(int));
      Q->max_position = (int *)malloc(Q->size * sizeof(int));
      Q->a = (double **)malloc(Q->size * sizeof(double *));
      for (int i = 0; i < Q->size; i++){
         Q->position[i] = -1;
         Q->max_position[i] = q_lens[i]-1;
         Q->a[i] = (double *)malloc(q_lens[i] * sizeof(double));
      }
   }
   else {
      Q->q = new queue<double>[Q->size];
   }
}

void qFree(Queue *Q)
{
   if (Q->type == Q_ARRAY){
      for (int i = 0; i < Q->size; i++){
         free(Q->a[i]);
      }
      free(Q->a);
      free(Q->position);
      free(Q->max_position);
   }
   else {
      delete [] Q->q;
   }
}

void qPut(Queue *Q, int destinationQID, double sourceData)
{
   int i = destinationQID;

   omp_set_lock(&(Q->lock));

   if (Q->type == Q_ARRAY){   
      if (Q->position[i] < Q->max_position[i]){
         Q->position[i]++;
         Q->a[i][Q->position[i]] = sourceData;
      }
   }
   else {
      Q->q[i].push(sourceData);
   }

   omp_unset_lock(&(Q->lock));
}

int qPoll(Queue *Q, int destinationQID, double *destinationData)
{
   int flag = 0;
   int i = destinationQID;
   double x;
    
   omp_set_lock(&(Q->lock));

   if (Q->type == Q_ARRAY){
      if (Q->position[i] > -1){
         *destinationData = Q->a[i][0];
         flag = 1;
      }
   }
   else {
      if (!(Q->q[i].empty())){
         *destinationData = Q->q[i].front();
         flag = 1;
      }
   }

   omp_unset_lock(&(Q->lock));

   return flag;
}

void qWait(Queue *Q, int destinationQID, double *destinationData)
{
   int q_empty = 1;
   int i = destinationQID;
   double x;

   while (qPoll(Q, i, &x) == 0);

   *destinationData = x;
}

int qGet(Queue *Q, int destinationQID, double *destinationData)
{
   int q_empty = 1;
   int i = destinationQID;
   double x;
   int flag = 0;

   omp_set_lock(&(Q->lock));

   if (Q->type == Q_ARRAY){
      if (Q->position[i] > -1){
         *destinationData = Q->a[i][0];
         qPopArray(Q, i);
         Q->position[i]--;
         flag = 1;
      }
   }
   else {
      if (!(Q->q[i].empty())){
         *destinationData = Q->q[i].front();
         Q->q[i].pop();
         flag = 1;
      }
   }

   omp_unset_lock(&(Q->lock));

   return flag;
}

void qAccum(Queue *Q, int destinationQID, double y)
{
   double x;
   qWait(Q, destinationQID, &x);
   x += y;
   qPut(Q, destinationQID, x);
   qGet(Q, destinationQID, &x);
}
