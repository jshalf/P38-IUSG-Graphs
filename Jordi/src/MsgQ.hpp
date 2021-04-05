#ifndef MSGQ_HPP
#define MSGQ_HPP

#include "Main.hpp"

/* queue data structure */
typedef struct {
   queue<double> *q; /* std queue */
   int size = 0; /* size of queue */
   omp_lock_t *lock; /* OpenMP lock */
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

Queue* qGetObj(void);

#endif
