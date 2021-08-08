#ifndef JACOBI_HPP
#define JACOBI_HPP

#include "../../src/Main.hpp"
#include "../../src/Matrix.hpp"

void Histogram_Seq(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   );

void Histogram_Par(HistogramData *hist,
                   int *index, int n_index,
                   int *Tally, int n_tally
                   );

#endif
