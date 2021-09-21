#ifndef MSGQ_HPP
#define MSGQ_HPP

#include "Main.hpp"

#define PROGRESS_QPUT 0
#define PROGRESS_QGET 1

typedef struct {
   queue<double> *wtime = NULL;
   queue<int> *idx = NULL;
   queue<double> *wtime_glob = NULL;
   double wtime_glob_start;
   int size;
} DAG;

/* queue data structure */
typedef struct {
   queue<double> *q; /* std queue */
   int size = 0; /* size of queue */
   omp_lock_t *lock; /* OpenMP lock */
   DAG dag;
   int dag_flag = 0;
} Queue;

void qInitLock(Queue *Q);

void qDestroyLock(Queue *Q);

void qAlloc(Queue *Q, int n);

void qFree(Queue *Q);

void qPut(Queue *Q, int destinationQID, double sourceData);

int qPoll(Queue *Q, int destinationQID, double *destinationData);

void qWait(Queue *Q, int destinationQID, double *destinationData);

int qGet(Queue *Q, int destinationQID, double *destinationData);

void qAccum(Queue *Q, int destinationQID, double y);

double qGetWtime(void);

Queue* qGetObj(void);

void qPrintDAG(Queue *Q);

#endif
