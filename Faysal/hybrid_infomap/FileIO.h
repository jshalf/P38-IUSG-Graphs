
#ifndef FILEIO_H
#define FILEIO_H

#include <iostream>
#include <sstream>
#include <fstream>
#include "Node.h"
#include "Module.h"
#include </global/common/sw/cray/cnl7/haswell/metis/5.1.0/intel/19.0.3.199/nbtsmmb/include/metis.h>

using namespace std;

// declaration of functions related to file IO.
void load_pajek_format_network(string fName, Network &network);
void load_linkList_format_network(string fName, Network &network);
void load_csr_format_network(string fName, Network &network, string metisFile);
void old_load_csr_format_network(string fName, Network &network);
void write_cluster_result(string outFile, Network &finalNetwork);

#endif
