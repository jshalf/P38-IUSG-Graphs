#ifndef MSGQ_HPP
#define MSGQ_HPP

#include <vector>
#include <queue>

#define PROGRESS_QPUT 0
#define PROGRESS_QGET 1

#define MSGQ_CACHE_LINE_SIZE 64

#ifdef USE_DEVA

#include "deva_includes.hpp"

using deva::process_n;
using deva::rank_n;
using deva::rank_me;
using deva::rank_me_local;
constexpr int proc_rank_n = rank_n / process_n;

#else
#include <omp.h>
#endif

template<typename T>
struct alignas(MSGQ_CACHE_LINE_SIZE) MsgQCacheAlign {
   T data;
};

// TODO: make a "queue" class that holds the message queue data and would use an std::queue underneath

template<typename MessageData>
class MessageQueue {
   protected:
      std::vector<std::queue<MessageData>> Qdata;
#ifdef USE_STDTHREAD
      std::mutex* lock;
#else
      std::vector<MsgQCacheAlign<omp_lock_t>> lock;
#endif
      void qInitLock(int n);
      void qDestroyLock(void);
      void qAlloc(int n);
      void qFree(void);

#if USE_DEVA
      std::queue<MessageData>& GetMyQ(void);
#endif

   public:
      MessageQueue(int n);
      ~MessageQueue(void);

      /* message queue "put" primitive */
      void qPut(int destinationQID,      /* destination queue id */
                MessageData sourceData); /* source data */
      /* message queue "get" primitive */
      int qGet(int destinationQID,            /* destination queue id */
               MessageData *destinationData); /* destination data */
      /* message queue "poll" primitive */
      int qPoll(int destinationQID,            /* destination queue id */
                MessageData *destinationData); /* destination data */
      
};

template<typename MessageData>
MessageQueue<MessageData>::MessageQueue(int n)
{
   qAlloc(n);
   qInitLock(n);
}

template<typename MessageData>
MessageQueue<MessageData>::~MessageQueue(void)
{
   qDestroyLock();
   qFree();
}

/* initialize mutex lock */
template<typename MessageData>
void MessageQueue<MessageData>::qInitLock(int n)
{
#ifdef USE_STDTHREAD
#elif USE_DEVA
#else
   for (int i = 0; i < lock.size(); i++){
      omp_init_lock(&(lock[i].data));
   }
#endif
}

/* destroy mutex  */
template<typename MessageData>
void MessageQueue<MessageData>::qDestroyLock(void)
{
#ifdef USE_STDTHREAD
#elif USE_DEVA
#else
   for (int i = 0; i < lock.size(); i++){
      omp_destroy_lock(&(lock[i].data));
   }
#endif
}

/* initialize queue data */
template<typename MessageData>
void MessageQueue<MessageData>::qAlloc(int n) /* length of vector of queues */
{
   Qdata.resize(n);
#ifdef USE_STDTHREAD
   lock = (std::mutex *)malloc(sizeof(std::mutex) * n);
#elif USE_DEVA
#else
   lock.resize(n);
#endif
}

/* free queue data */
template<typename MessageData>
void MessageQueue<MessageData>::qFree(void)
{
#ifdef USE_STDTHREAD
   free(lock);
#elif USE_DEVA
#else
#endif
}

#if USE_DEVA
template<typename MessageData>
std::queue<MessageData>& MessageQueue<MessageData>::GetMyQ(void)
{
  return Qdata[rank_me_local()];
}
#endif

/*****************
 * put primitive 
 *****************/
template<typename MessageData>
void MessageQueue<MessageData>::qPut(int destinationQID,       /* destination queue id */
                                     MessageData sourceData)  /* source data */
{
/* If using pthreads or OpenMP, acquire lock.
 * If using devastator, send lambda start here. */
#ifdef USE_STDTHREAD
   lock[destinationQID].data();
#elif USE_DEVA
   deva::send(destinationQID,
      [this] (MessageData Msg) { 
#else
   omp_set_lock(&(lock[destinationQID].data)); /* acquire lock */
#endif




/* push to queue */
#ifdef USE_DEVA
         //Qdata[destinationQID].push(std::move(Msg));
         GetMyQ().push(std::move(Msg));
#else
   Qdata[destinationQID].push(sourceData);
#endif





/* If using pthreads or OpenMP, release lock.
 * If using devastator, lambda ends here. */
#ifdef USE_STDTHREAD
   lock[destinationQID].unlock();
#elif USE_DEVA
      }, 
      std::move(sourceData)
   );
#else
   omp_unset_lock(&(lock[destinationQID].data));
#endif
}

/* message queue poll primitive */
template<typename MessageData>
int MessageQueue<MessageData>::qPoll(int destinationQID,           /* destination queue id */
                                     MessageData *destinationData) /* destination data */
{
   int flag = 0;
   double x;
#ifdef USE_STDTHREAD
   lock[destinationQID].data();
#else
   omp_test_lock(&(lock[destinationQID].data));
#endif

   /* check for messages in the queue */
   if (Qdata[destinationQID].size()){
      /* if queue is not empty, get front and flag of one */
      *destinationData = Qdata[destinationQID].front();
      flag = 1;
   }

#ifdef USE_STDTHREAD
   lock[destinationQID].unlock();
#else
   omp_unset_lock(&(lock[destinationQID].data)); /* release lock */
#endif

   return flag;
}

/* message queue get primitive */
template<typename MessageData>
int MessageQueue<MessageData>::qGet(int destinationQID,           /* destination queue id */
                                    MessageData *destinationData) /* destination data */
{
   double x;
   int flag = 0;
/* acquire lock */
#ifdef USE_STDTHREAD
   lock[destinationQID].data();
#elif USE_DEVA
#else
   omp_set_lock(&(lock[destinationQID].data));
#endif

   /* Check for messages in the queue */
   if (Qdata[destinationQID].size()){
      /* if queue is not empty, pop front and flag of one */
#ifdef USE_DEVA
      *destinationData = std::move(Qdata[destinationQID].front());
#else
      *destinationData = Qdata[destinationQID].front();
#endif
      Qdata[destinationQID].pop(); /* pop front */
      flag = 1;
   }

/* release lock */
#ifdef USE_STDTHREAD
   lock[destinationQID].unlock();
#elif USE_DEVA
#else
   omp_unset_lock(&(lock[destinationQID].data));
#endif

   return flag;
}

#endif
