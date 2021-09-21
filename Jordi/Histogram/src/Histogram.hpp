#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

typedef struct {
   int thread_recv_count;
   int actual_num_qGets;
   int num_qGets;
   int *num_qPuts;
   int tid;
   int num_threads;
} HistogramProgressData;

void Histogram_Seq(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   );

void Histogram_Par(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   );

#endif
