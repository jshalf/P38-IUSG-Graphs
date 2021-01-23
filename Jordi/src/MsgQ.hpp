#ifndef MSGQ_HPP
#define MSGQ_HPP

#define Q_ARRAY 0
#define Q_STDQUEUE 1

#include "Main.hpp"

typedef struct {
   double **a;
   queue<double> *q;
   int *max_position;
   int *position;
   int size;
   int type;
   omp_lock_t lock;
} Queue;

void qInitLock(Queue *Q);

void qDestroyLock(Queue *Q);

//void qAlloc(int n);
void qAlloc(Queue *Q, int n, int *q_lens);

void qFree(Queue *Q);

void qPut(Queue *Q, int destinationQID, double sourceData);

int qPoll(Queue *Q, int destinationQID, double *destinationData);

void qWait(Queue *Q, int destinationQID, double *destinationData);

int qGet(Queue *Q, int destinationQID, double *destinationData);

void qAccum(Queue *Q, int destinationQID, double y);

Queue* qGetObj(void);
//queue<double>* qGetObj(void);

#endif
