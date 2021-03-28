#include <iostream>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include "Module.h"
#include "Node.h"
#include "timing.h"
#include <mpi.h>
#include <cstdio>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <malloc.h>
#include <string>
#include <sstream>
#include "global.h"
#include <unordered_map>
#include </opt/intel/vtune_profiler_2020.2.0.610396/include/ittnotify.h>
#include <chrono>


//typedef __gnu_cxx::hash_map<int, double> flowmap;
typedef unordered_map<int, double> flowmap;
typedef map<int, pair<double, double> > modInfo;	// <modID, <exitPr, sumPr> >

void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);

const double InitLength = 10000.0;

using namespace std;

struct MoveSummary {
	double diffCodeLen;			// delta(L(M))
	int newModule;
	double sumPr1, sumPr2;// updated sumPr1 and sumPr2.  1 --> oldM, 2 --> newM.
	double exitPr1, exitPr2;	// updated exitPr1 and exitPr2.
	double newSumExitPr;// SUM(q_i) = this.sumAllExitPr + q_1_new + q_2_new - q_1_old - q_2_old.
};

const double Network::alpha;
const double Network::beta;
// Constructors of Module class.
Module::Module() :
		index(0), exitPr(0.0), stayPr(0.0), sumPr(0.0), sumTPWeight(0.0), sumDangling(
				0.0), numMembers(1) {
}

Module::Module(int idx, double exitPr, double sumPr) :
		index(idx), exitPr(exitPr), stayPr(exitPr + sumPr), sumPr(sumPr), sumTPWeight(
				0.0), sumDangling(0.0), numMembers(1) {
}

//Module::Module(int idx, int nNode, Node nd)
Module::Module(int idx, Node * nd) :
		index(idx), exitPr(0.0), stayPr(0.0), sumPr(nd->Size()), sumTPWeight(
				nd->TeleportWeight()), sumDangling(0.0), numMembers(1) {
	double sumExitFlow = 0.0;

	int nOutLinks = nd->outLinks.size();

	for (int i = 0; i < nOutLinks; i++)
		sumExitFlow += nd->outLinks[i].second;// after call initiate(), w_ab is updated as p_a * w_ab.

	if (nd->IsSuper()) {// If SuperNode, we already calculated exitPr in the corresponding module.
		exitPr = nd->ExitPr();
		sumDangling = nd->DanglingSize();
	}
	// exitPr = tau * (1 - sum(tau_a))*sum(p_a) + (1-tau) * sumExitFlow.
	else if (!nd->IsDangling()) {
		exitPr = Network::alpha * (1.0 - nd->TeleportWeight()) * sumPr
				+ Network::beta * sumExitFlow;
		nd->setExitPr(exitPr);
	} else {
		exitPr = (1.0 - nd->TeleportWeight()) * sumPr;
		nd->setExitPr(exitPr);
		sumDangling = nd->Size();
	}

	stayPr = exitPr + sumPr;	// exitPr + sum (P_a).
	members.push_back(nd);
}

SubModule::SubModule() :
		modIdx(0), numMembers(0), sumPr(0.0), sumTPWeight(0.0), sumDangling(0.0) {
}

/**
 *	Module mod : module of newNetwork (network based on a module in original network.)
 */
SubModule::SubModule(Module& mod, map<int, int>& origNodeID, int modIndex) {
	numMembers = mod.numMembers;
	sumPr = mod.sumPr;
	sumTPWeight = mod.sumTPWeight;
	sumDangling = mod.sumDangling;
	modIdx = modIndex;

	for (vector<Node*>::iterator it = mod.members.begin();
			it != mod.members.end(); it++) {
		members.push_back(origNodeID[(*it)->ID()]);
	}

	// numMembers should be equal to members.size() ...
	if (numMembers != members.size())
		cout
				<< "SOMETHING WRONG!!! -- numMembers != members.size() in SubModule()..."
				<< endl;
}

SubModule::SubModule(Module& mod) {
	modIdx = mod.index;
	numMembers = 1;
	sumPr = mod.sumPr;
	sumTPWeight = mod.sumTPWeight;
	sumDangling = mod.sumDangling;

	members.push_back(mod.members[0]->ID());
}

// Constructors of Network class.
Network::Network() :
		level(0), codeLength(InitLength), nNode(0), nEdge(0), nEmptyMod(0), nModule(
				0), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
}

Network::Network(int numNode, int numEdge, int level) :
		level(level), codeLength(InitLength), nNode(numNode), nEdge(numEdge), nEmptyMod(
				0), nModule(numNode), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
	emptyModules.reserve(numNode);
}

Network::Network(int numNode, int numEdge, int level, double codeLen) :
		level(level), codeLength(codeLen), nNode(numNode), nEdge(numEdge), nEmptyMod(
				0), nModule(numNode), totNodeWeights(0.0), nDanglings(0), allNodes_log_allNodes(
				0.0), sumAllExitPr(0.0) {
	emptyModules.reserve(numNode);
}

/*
 * Other Member functions of Network class.
 */

void Network::findDanglingNodes() 
{

	nDanglings = 0;		// reset the number of dangling nodes.

	for (int i = 0; i < nNode; i++) {
		if (nodes[i].outLinks.empty()) {
			danglings.push_back(i);
			nDanglings++;
			nodes[i].setIsDangling(true);
		}
	}

}

/*
//void Network::initiate() {
void Network::initiate(int numTh) 
{
	int size, rank;

	struct timeval startIni, endIni;

	omp_set_num_threads(numTh);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	const int processor_counts = size;

	int nDangNodes = 0;

	struct timeval startT, endT;
	
	auto tik = std::chrono::high_resolution_clock::now();

	gettimeofday(&startT, NULL);

	#pragma omp parallel reduction (+:nDangNodes)
	 {
	 int myID = omp_get_thread_num();
	 int nTh = omp_get_num_threads();

	 int startIndex, endIndex;
	 findAssignedPart(&startIndex, &endIndex, nNode, nTh, myID);

	for (int i = 0; i < nNode; i++) 
	{
		if (nodes[i].outLinks.empty()) 
		{
			#pragma omp critical (updateDang)
			 {
			danglings.push_back(i);
			//}
			nDangNodes++;
			nodes[i].setIsDangling(true);
		} 
		else 
		{	// normal nodes. --> Normalize edge weights.
			int nOutLinks = nodes[i].outLinks.size();
			//double sum = nodes[i].selfLink;	// don't support selfLink yet.
			double sum = 0.0;
			for (int j = 0; j < nOutLinks; j++) 
			{
				sum += nodes[i].outLinks[j].second;
			}
			for (int j = 0; j < nOutLinks; j++) 
			{
				nodes[i].outLinks[j].second /= sum;
			}
		}
	}
	//}

	gettimeofday(&endT, NULL);

	nDanglings = nDangNodes;

	cout << "Level " << level << ": the number of dangling nodes = "
			<< nDanglings << endl;
	cout << "Time for finding dangling nodes : "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);

	calculateSteadyState(numTh);

	gettimeofday(&endT, NULL);
	cout << "Time for calculating steady state of nodes (eigenvector): "
			<< elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);

	// Update edges to represent flow based on calculated steady state (aka size).
	for (int i = 0; i < nNode; i++) 
	{
		if (!nodes[i].IsDangling()) 
		{
			int nOutLinks = nodes[i].outLinks.size();
			for (int j = 0; j < nOutLinks; j++) 
			{
				nodes[i].outLinks[j].second = nodes[i].Size() * nodes[i].outLinks[j].second;
			}
		} 
		else 
		{
			nodes[i].setDanglingSize(nodes[i].Size());
		}
	}

	gettimeofday(&endT, NULL);
	cout << "Time for updating edge weights based on flow information: " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	// Calculate SUM of p_log_p over all nodes.
	int tag = 888888;
	double allNds_log_allNds_send = 0.0;
	double allNds_log_allNds_recv = 0.0;
	double* allNds_log_allNds_recv_all = new double[processor_counts]();

	int start, end;
	findAssignedPart(&start, &end, nNode, size, rank);

	for (int i = start; i < end; i++) 
	{
		allNds_log_allNds_send += pLogP(nodes[i].Size());
	}

	allNodes_log_allNodes = allNds_log_allNds_send;

	gettimeofday(&startIni, NULL);
	MPI_Allgather(&allNds_log_allNds_send, 1, MPI_DOUBLE, allNds_log_allNds_recv_all, 1, MPI_DOUBLE, MPI_COMM_WORLD);

	for (int p = 0; p < size; p++) 
	{
		if (p != rank) 
		{
			allNodes_log_allNodes += allNds_log_allNds_recv_all[p];
		}
	}

	gettimeofday(&endIni, NULL);

	/////////////////////////////
	// Make modules from nodes //
	/////////////////////////////
	for (int i = 0; i < nNode; i++) 
	{
		modules[i] = Module(i, &nodes[i]);// Assign each Node to a corresponding Module object. Initially, each node is in its own module.
		nodes[i].setModIdx(i);
	}

	nModule = nNode;

	calibrate(numTh, 0); //sending tag 0 to indicate the calibrate call from initiate function

	delete[] allNds_log_allNds_recv_all;

	auto tok = std::chrono::high_resolution_clock::now();

	initiationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

}
*/

void Network::initiate(int numTh) 
{

	int nDangNodes = 0;

	auto tik = std::chrono::high_resolution_clock::now();

	struct timeval startT, endT;

	omp_set_num_threads(numTh);

	gettimeofday(&startT, NULL);

	#pragma omp parallel reduction (+:nDangNodes) 
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) 
		{
			if (nodes[i].outLinks.empty()) 
			{
				#pragma omp critical (updateDang)
				{
					danglings.push_back(i);
				}
				nDangNodes++;
				nodes[i].setIsDangling(true);
			}
			else 
			{	// normal nodes. --> Normalize edge weights.
				int nOutLinks = nodes[i].outLinks.size();
				//double sum = nodes[i].selfLink;	// don't support selfLink yet.
				double sum = 0.0;
				for (int j = 0; j < nOutLinks; j++)
				{
					sum += nodes[i].outLinks[j].second;
				}
				//nodes[i].selfLink /= sum;
				for (int j = 0; j < nOutLinks; j++)
				{
					nodes[i].outLinks[j].second /= sum;
				}
			}
		}
	}

	gettimeofday(&endT, NULL);

	nDanglings = nDangNodes;

	cout << "Level " << level << ": the number of dangling nodes = " << nDanglings << endl;
	cout << "Time for finding dangling nodes : " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;

	gettimeofday(&startT, NULL);
	calculateSteadyState(numTh);
	gettimeofday(&endT, NULL);
	cout << "Time for calculating steady state of nodes (eigenvector): " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;


	gettimeofday(&startT, NULL);

	// Update edges to represent flow based on calculated steady state (aka size).
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);

		for (int i = start; i < end; i++) 
		{
			if (!nodes[i].IsDangling()) 
			{
				int nOutLinks = nodes[i].outLinks.size();
				for (int j = 0; j < nOutLinks; j++) 
				{
					nodes[i].outLinks[j].second = nodes[i].Size() * nodes[i].outLinks[j].second;
				}
			}
			else 
			{
				nodes[i].setDanglingSize(nodes[i].Size());
			}
		}
	}

	gettimeofday(&endT, NULL);
	cout << "Time for updating edge weights based on flow information: " << elapsedTimeInSec(startT, endT) << " (sec)" << endl;
	
	// Calculate SUM of p_log_p over all nodes.
	double allNds_log_allNds = 0.0;

	#pragma omp parallel reduction (+:allNds_log_allNds)
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);
	
		for (int i = start; i < end; i++)
		{
			allNds_log_allNds += pLogP(nodes[i].Size());
		}
	}

	allNodes_log_allNodes = allNds_log_allNds;

	/////////////////////////////
	// Make modules from nodes //
	/////////////////////////////
	#pragma omp parallel
	{
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;
		findAssignedPart(&start, &end, nNode, nTh, myID);
		for (int i = start; i < end; i++) 
		{
			modules[i] = Module(i, &nodes[i]);	// Assign each Node to a corresponding Module object.
			// Initially, each node is in its own module.
			nodes[i].setModIdx(i);
		}
	}

	nModule = nNode;
	calibrate(numTh, 0);
	auto tok = std::chrono::high_resolution_clock::now();
	initiationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();
}


// calculating steady state of nodes by Power Iteration Method.
// same as Greedy::eigenvector() in Infomap implementation.
// Modify for the design of this implementation.

void Network::calculateSteadyState(int numTh) 
{
	
	// initial probability distribution = 1 / N.
	vector<double> size_tmp = vector<double>(nNode, 1.0 / nNode);

	int iter = 0;
	double danglingSize = 0.0;
	double sqdiff = 1.0;
	double sum = 0.0;
	
	auto tik = std::chrono::high_resolution_clock::now();

	// Generate addedSize array per each thread, so that we don't need to use lock for each nodeSize addition.
	double** addedSize = new double*[numTh];
	omp_set_num_threads(numTh);

	for (int i = 0; i < numTh; i++)
	{
		addedSize[i] = new double[nNode]();
	}

	do {
		// calculate sum of the size of dangling nodes.
		danglingSize = 0.0;

#pragma omp parallel reduction (+:danglingSize)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nDanglings, nTh, myID);

			for (int i = start; i < end; i++)
			{
				danglingSize += size_tmp[danglings[i]];
			}
		}

		// flow via teleportation.
#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++)
			{
				nodes[i].setSize((alpha + beta * danglingSize)* nodes[i].TeleportWeight());
			}
		}

		int realNumTh = 0;

		// flow from network steps via following edges.
#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

#pragma omp master
			{
				realNumTh = nTh;
			}

			double *myAddedSize = addedSize[myID];

			for (int i = 0; i < nNode; i++)
			{
				myAddedSize[i] = 0.0; // initialize for the temporary addedSize array.
			}

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++) 
			{
				int nOutLinks = nodes[i].outLinks.size();

				for (int j = 0; j < nOutLinks; j++) 
				{
					myAddedSize[nodes[i].outLinks[j].first] += beta * nodes[i].outLinks[j].second * size_tmp[i];
				}
			}
		}
		//set nodes[i].size by added addedSize[] values.
#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++) 
			{
				for (int j = 0; j < realNumTh; j++)
				{
					nodes[i].addSize(addedSize[j][i]);
				}
			}
		}

		// Normalize of node size.
		sum = 0.0;
#pragma omp parallel reduction (+:sum)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++)
			{
				sum += nodes[i].Size();
			}
		}
		sqdiff = 0.0;

#pragma omp parallel reduction (+:sqdiff)
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++) 
			{
				nodes[i].setSize(nodes[i].Size() / sum);
				sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
				size_tmp[i] = nodes[i].Size();
			}
		}

		iter++;

	} while ((iter < 200) && (sqdiff > 1.0e-15 || iter < 50));

	// deallocate 2D array.
	for (int i = 0; i < numTh; i++)
	{
		delete[] addedSize[i];
	}
	delete[] addedSize;

	auto tok = std::chrono::high_resolution_clock::now();

	computePageRankTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();
}


// This is pagerank calculation with no explicit work division by multithreading
// Rather the omp itself takes care of the work load division

/*
void Network::calculateSteadyState(int numTh) {
	
	// initial probability distribution = 1 / N.

	vector<double> size_tmp = vector<double>(nNode, 1.0 / nNode);

	int iter = 0;
	double danglingSize = 0.0;
	double sqdiff = 1.0;
	double sum = 0.0;
	
	auto tik = std::chrono::high_resolution_clock::now();

	// Generate addedSize array per each thread, so that we don't need to use lock for each nodeSize addition.
	double** addedSize = new double*[numTh];
	omp_set_num_threads(numTh);

	for (int i = 0; i < numTh; i++)
	{
		addedSize[i] = new double[nNode]();
	}

	do {
		// calculate sum of the size of dangling nodes.
		danglingSize = 0.0;

#pragma omp parallel for reduction (+:danglingSize)

			for (int i = 0; i < nDanglings; i++)
			{
				danglingSize += size_tmp[danglings[i]];
			}


		// flow via teleportation.
#pragma omp parallel for
			for (int i = 0; i < nNode; i++)
			{
				nodes[i].setSize((alpha + beta * danglingSize)* nodes[i].TeleportWeight());
			}

		int realNumTh = 0;

		// flow from network steps via following edges.
#pragma omp parallel
		{
			int myID = omp_get_thread_num();
			int nTh = omp_get_num_threads();

#pragma omp master
			{
				realNumTh = nTh;
			}

			double *myAddedSize = addedSize[myID];

			for (int i = 0; i < nNode; i++)
			{
				myAddedSize[i] = 0.0; // initialize for the temporary addedSize array.
			}

			int start, end;
			findAssignedPart(&start, &end, nNode, nTh, myID);

			for (int i = start; i < end; i++) 
			{
				int nOutLinks = nodes[i].outLinks.size();

				for (int j = 0; j < nOutLinks; j++) 
				{
					myAddedSize[nodes[i].outLinks[j].first] += beta * nodes[i].outLinks[j].second * size_tmp[i];
				}
			}
		}

		//set nodes[i].size by added addedSize[] values.
#pragma omp parallel for
			for (int i = 0; i < nNode; i++) 
			{
				for (int j = 0; j < realNumTh; j++)
				{
					nodes[i].addSize(addedSize[j][i]);
				}
			}

		// Normalize of node size.
		sum = 0.0;
#pragma omp parallel for reduction (+:sum)
			for (int i = 0; i < nNode; i++)
			{
				sum += nodes[i].Size();
			}

		sqdiff = 0.0;

#pragma omp parallel for reduction (+:sqdiff)
			for (int i = 0; i < nNode; i++) 
			{
				nodes[i].setSize(nodes[i].Size() / sum);
				sqdiff += fabs(nodes[i].Size() - size_tmp[i]);
				size_tmp[i] = nodes[i].Size();
			}

		iter++;

	} while ((iter < 200) && (sqdiff > 1.0e-15 || iter < 50));

	// deallocate 2D array.
	for (int i = 0; i < numTh; i++)
	{
		delete[] addedSize[i];
	}
	delete[] addedSize;


	auto tok = std::chrono::high_resolution_clock::now();

	computePageRankTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

	cout << "Calculating flow done in " << iter << " iterations!" << endl;
}
*/

// This function calculate current codeLength.
// This implementation is modified version of infomap implementation.
void Network::calibrate(int numTh, int tag) {
	//This is the calculation of Equation (4) in the paper.

	int size;
	int rank;

	int arraySize = 3;

	auto tik = std::chrono::high_resolution_clock::now();
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double* send_exit_stay_sum = new (nothrow) double[arraySize]();
	//double receive_sum_exit_stay[arraySize] = { 0.0, 0.0, 0.0 };
	//double* s_sum_exit_stay = new (nothrow) double[arraySize]();
	double* receive_exit_stay_sum = new (nothrow) double[arraySize * size]();

	if (!receive_exit_stay_sum || !send_exit_stay_sum) {
		printf(
				"ERROR in calibrate: process:%d does not have sufficient memory, process exiting...\n",
				rank);
		exit(1);
	}

	/*	double send_exit_stay_sum[3] = { 0.0, 0.0, 0.0 };
	 double receive_exit_stay_sum[3] = { 0.0, 0.0, 0.0 };*/

	double sum_exit_log_exit = 0.0;
	double sum_stay_log_stay = 0.0;
	double sumExit = 0.0;

	int start, end;
	findAssignedPart(&start, &end, nModule, size, rank);

	for (int i = start; i < end; i++) {
		sum_exit_log_exit += pLogP(modules[i].exitPr);
		sum_stay_log_stay += pLogP(modules[i].stayPr);
		sumExit += modules[i].exitPr;
	}

	send_exit_stay_sum[0] = sum_exit_log_exit;
	send_exit_stay_sum[1] = sum_stay_log_stay;
	send_exit_stay_sum[2] = sumExit;

	MPI_Allgather(send_exit_stay_sum, arraySize, MPI_DOUBLE,
			receive_exit_stay_sum, arraySize, MPI_DOUBLE,
			MPI_COMM_WORLD);

	for (int p = 0; p < size; p++) {
		if (rank != p) {
			int startIndex = p * arraySize;
			sum_exit_log_exit += receive_exit_stay_sum[startIndex];
			sum_stay_log_stay += receive_exit_stay_sum[startIndex + 1];
			sumExit += receive_exit_stay_sum[startIndex + 2];
		}
	}

	sumAllExitPr = sumExit;
	double sumExit_log_sumExit = pLogP(sumExit);

	codeLength = sumExit_log_sumExit - 2.0 * sum_exit_log_exit
			+ sum_stay_log_stay - allNodes_log_allNodes;
	
	auto tok = std::chrono::high_resolution_clock::now();

	calibrationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

	delete[] send_exit_stay_sum;
	delete[] receive_exit_stay_sum;
}

/*
 * This function implements the 1) procedure of stochastic_greedy_partition(Network &network) function in SeqInfomap.cpp.
 *
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 *	
 *	return: the number of moves.
 */
int Network::move(int iteration) 
{

	int size;
	int rank;
	int elementCount = 0;
	int mpi_tag = 3333 + iteration;
	srand(time(0));

	struct timeval startMove, endMove;

	gettimeofday(&startMove, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// Generate random sequential order of nodes.

	double totalElements = nNode;

	int stripSize = ceil(totalElements / size) + 1;

	const unsigned int intPackSize = 2 * stripSize + 3; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 4 * stripSize + 3;

	int* randomGlobalArray = new int[nNode]();

	int start, end;
	findAssignedPart(&start, &end, nNode, size, rank);

	int numberOfElements = end - start; //number of elements assigned to this processor

	for (int i = start; i < end; i++) 
	{
		randomGlobalArray[i] = i;
	}

	for (int i = start; i < end; i++) 
	{
		int pretarget = rand() % numberOfElements;
		int target = start + pretarget;
		// swap numbers between i and target.
		int tmp = randomGlobalArray[i];
		randomGlobalArray[i] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}

	int numMoved = 0;	// the counter for the number of movements.

	int* intSendPack = new (nothrow) int[intPackSize]();// multiplier 3 for 3 different elements of index, newModules, oldModules and +1 for elementCount
	int* intReceivePack = new (nothrow) int[intPackSize]();
	int* intAllRecvPack = new (nothrow) int[intPackSize * size]();
	double* doubleSendPack = new (nothrow) double[doublePackSize]();// multiplier 6 for 6 different elements of diffCodeLen, sumPr1, sumPr2, exitPr1, exitPr2, newSumExitPr
	double* doubleReceivePack = new (nothrow) double[doublePackSize]();
	double* doubleAllRecvPack = new (nothrow) double[doublePackSize * size]();

	if (!intSendPack || !intReceivePack || !intAllRecvPack || !doubleSendPack
			|| !doubleReceivePack || !doubleAllRecvPack) {
		printf(
				"ERROR in move: process:%d does not have sufficient memory, process exiting...\n",
				rank);
		exit(1);
	}

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	elementCount = 0;

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = start; i < end; i++) 
	{

		int counter = 0;
		bool moduleUpdate = false;

		int vertexIndex = randomGlobalArray[i];
		Node& nd = nodes[vertexIndex];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		// During this counting, aggregate outFlowNotToOldMod and inFlowFromOldMod values.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();		// p_nd.
		double ndTPWeight = nd.TeleportWeight();		// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		int newMod;
		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod;

		counter++;

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.

					// update the moduleUpdate flag to indicate that an update has been made for the current vertex i
					moduleUpdate = true;

					//save new module id
					bestResult.newModule = currentResult.newModule;

					bestResult.diffCodeLen = currentResult.diffCodeLen;

					bestResult.sumPr1 = currentResult.sumPr1;

					bestResult.sumPr2 = currentResult.sumPr2;

					bestResult.exitPr1 = currentResult.exitPr1;

					bestResult.exitPr2 = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			}
		}

		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
			// update related to newMod...
			int newMod = bestResult.newModule;

			if (oldMod < newMod) {

				if (modules[newMod].numMembers == 0) {
					nEmptyMod--;
					nModule++;
				}

				nd.setModIdx(newMod);

				modules[newMod].numMembers++;
				modules[newMod].exitPr = bestResult.exitPr2;
				modules[newMod].sumPr = bestResult.sumPr2;
				modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
				modules[newMod].sumTPWeight += ndTPWeight;

				if (nd.IsDangling()) {
					modules[newMod].sumDangling += ndSize;
					modules[oldMod].sumDangling -= ndSize;
				}

				// update related to the oldMod...
				modules[oldMod].numMembers--;
				modules[oldMod].exitPr = bestResult.exitPr1;
				modules[oldMod].sumPr = bestResult.sumPr1;
				modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
				modules[oldMod].sumTPWeight -= ndTPWeight;

				if (modules[oldMod].numMembers == 0) {
					nEmptyMod++;
					nModule--;
				}

				// the following code block is for sending information of updated vertex across the other processes

				intSendPack[numMoved] = vertexIndex;
				intSendPack[1 * stripSize + numMoved] = bestResult.newModule;

				doubleSendPack[numMoved] = bestResult.sumPr1;
				doubleSendPack[1 * stripSize + numMoved] = bestResult.sumPr2;
				doubleSendPack[2 * stripSize + numMoved] = bestResult.exitPr1;
				doubleSendPack[3 * stripSize + numMoved] = bestResult.exitPr2;

				sumAllExitPr = bestResult.newSumExitPr;

				codeLengthReduction += bestResult.diffCodeLen;

				numMoved++;
			}
		}
		elementCount++;
	}

	if (numberOfElements != elementCount) {
		printf(
				"something is not right, for rank:%d, numberOfElements:%d does not match elementCount:%d\n",
				rank, numberOfElements, elementCount);

	}

	intSendPack[intPackSize - 1] = numMoved;
	doubleSendPack[doublePackSize - 1] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[doublePackSize - 2] = codeLengthReduction;

	MPI_Allgather(intSendPack, intPackSize, MPI_INT, intAllRecvPack,
			intPackSize, MPI_INT, MPI_COMM_WORLD);

	MPI_Allgather(doubleSendPack, doublePackSize, MPI_DOUBLE, doubleAllRecvPack,
			doublePackSize, MPI_DOUBLE, MPI_COMM_WORLD);

	for (int p = 0; p < size; p++) {

		if (p != rank) {
			int startIntIndex = p * intPackSize;
			int endIntIndex = ((p + 1) * intPackSize) - 1;
			int totalElementSentFromSender = intAllRecvPack[endIntIndex];
			int startDblIndex = p * doublePackSize;
			int endDblIndex = (p + 1) * doublePackSize;
			sumAllExitPr += doubleAllRecvPack[endDblIndex - 1];
			codeLength += doubleAllRecvPack[endDblIndex - 2];
			int highest_index = totalElementSentFromSender + startIntIndex;

			for (int z = startIntIndex, d = startDblIndex; z < highest_index;
					z++, d++) {

				int nodeIndex = intAllRecvPack[z];

				Node& nod = nodes[nodeIndex];
				int oldModule = nod.modIdx;
				int newModule = intAllRecvPack[1 * stripSize + z];

				//we need to do some tweaking here to get rid of the extra reduction of codelength because of the circular moves in distributed informap
				//int nodeIndex = nod.ID();
				nod.setModIdx(newModule);

				if (oldModule != newModule) {

					numMoved++;

					if (modules[newModule].numMembers == 0) {
						nEmptyMod--;
						nModule++;
					}

					modules[newModule].numMembers++;
					modules[newModule].exitPr = doubleAllRecvPack[3 * stripSize
							+ d];
					modules[newModule].sumPr = doubleAllRecvPack[1 * stripSize
							+ d];
					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;

					modules[newModule].sumTPWeight += nod.TeleportWeight();

					if (nod.IsDangling()) {
						modules[newModule].sumDangling += nod.Size();
						modules[oldModule].sumDangling -= nod.Size();
					}

					modules[oldModule].numMembers--;
					modules[oldModule].exitPr = doubleAllRecvPack[2 * stripSize
							+ d];
					modules[oldModule].sumPr = doubleAllRecvPack[d];
					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;

					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					if (modules[oldModule].numMembers == 0) {
						nEmptyMod++;
						nModule--;
					}
				}
			}
		}
	}

	/*for (int processId = 0; processId < size; processId++) {
	 if (processId != rank) {
	 MPI_Sendrecv(intSendPack, intPackSize, MPI_INT, processId, mpi_tag,
	 intReceivePack, intPackSize, MPI_INT, processId, mpi_tag,
	 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

	 MPI_Sendrecv(doubleSendPack, doublePackSize, MPI_DOUBLE, processId,
	 mpi_tag, doubleReceivePack, doublePackSize, MPI_DOUBLE,
	 processId, mpi_tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

	 int totalElementSentFromSender = intReceivePack[intPackSize - 1];

	 sumAllExitPr += doubleReceivePack[doublePackSize - 1];
	 codeLength += doubleReceivePack[doublePackSize - 2];

	 for (int z = 0; z < totalElementSentFromSender; z++) {

	 int nodeIndex = intReceivePack[z];

	 Node& nod = nodes[nodeIndex];
	 int oldModule = nod.modIdx;
	 int newModule = intReceivePack[1 * stripSize + z];

	 //we need to do some tweaking here to get rid of the extra reduction of codelength because of the circular moves in distributed informap
	 //int nodeIndex = nod.ID();
	 nod.setModIdx(newModule);

	 if (oldModule != newModule) {
	 numMoved++;

	 if (modules[newModule].numMembers == 0) {
	 nEmptyMod--;
	 nModule++;
	 }

	 modules[newModule].numMembers++;
	 modules[newModule].exitPr = doubleReceivePack[3 * stripSize
	 + z];
	 modules[newModule].sumPr = doubleReceivePack[1 * stripSize
	 + z];
	 modules[newModule].stayPr = modules[newModule].exitPr
	 + modules[newModule].sumPr;

	 modules[newModule].sumTPWeight += nod.TeleportWeight();

	 if (nod.IsDangling()) {
	 modules[newModule].sumDangling += nod.Size();
	 modules[oldModule].sumDangling -= nod.Size();
	 }

	 modules[oldModule].numMembers--;
	 modules[oldModule].exitPr = doubleReceivePack[2 * stripSize
	 + z];
	 modules[oldModule].sumPr = doubleReceivePack[z];
	 modules[oldModule].stayPr = modules[oldModule].exitPr
	 + modules[oldModule].sumPr;
	 modules[oldModule].sumTPWeight -= nod.TeleportWeight();

	 if (modules[oldModule].numMembers == 0) {
	 nEmptyMod++;
	 nModule--;
	 }
	 }
	 }
	 }
	 }
	 */
	delete[] randomGlobalArray;
	delete[] intSendPack;
	delete[] intReceivePack;
	delete[] intAllRecvPack;
	delete[] doubleSendPack;
	delete[] doubleReceivePack;
	delete[] doubleAllRecvPack;

	gettimeofday(&endMove, NULL);

	return numMoved;
}

int Network::prioritize_move(double vThresh, int numTh, bool inWhile) 
{

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto start_time = std::chrono::high_resolution_clock::now();

	iteration_prioritizeMove++;

	const int totalProcessor = size;
	const int currentProcessor = rank;

	//mark 1 begin

	int *counts = new int[totalProcessor]();

	int *sendableElementCounts = new int[totalProcessor]();

	int *intCounts = new int[totalProcessor]();
	int *intDispls = new int[totalProcessor]();

	int *doubleCounts = new int[totalProcessor]();
	int *doubleDispls = new int[totalProcessor]();

	vector<int>().swap(myVertices); //initializing myVertices with an empty vector

	const int nActive = activeNodes.size();

	for (int i = 0; i < nActive; i++) {
		int originalNodeIndex = activeNodes[i];
		if (activeVertices.count(originalNodeIndex) > 0 && activeVertices.at(originalNodeIndex) == 1) 
		{
			counts[allNodes[originalNodeIndex]]++;

			if (allNodes[originalNodeIndex] == rank) 
			{
				myVertices.push_back(originalNodeIndex);
			}
		}
	}

	int myPortionVertexSize = myVertices.size(); //network.myVertices.size() basically equals to counts[rank]

	if (counts[rank] != myPortionVertexSize) 
	{
		printf("something is not right, counts[%d]:%d does not equal myPortionVertexSize:%d\n", rank, counts[rank], myPortionVertexSize);
	}
	const int intPackSize = 2 * myPortionVertexSize; // 1 extra space for storing the actual number of vertices that have been moved from one module to another
	const int doublePackSize = 2 * myPortionVertexSize; // 2 extra space, one for storing the sumAllExitPr, another for storing codeLengthReduction
	const int intAllPackSize = 2 * nActive;
	const int doubleAllPackSize = 2 * nActive;

	int* intSendPack = new (nothrow) int[intPackSize]();
	int* intAllRecvPack = new (nothrow) int[intAllPackSize](); // the extra size amount of memory is required for storing all the extra 1 space as 												  		//mentioned in intPacksize above

	double* doubleSendPack = new (nothrow) double[doublePackSize]();
	double* doubleAllRecvPack = new (nothrow) double[doubleAllPackSize](); // similar reason as above except the fact that we need extra 2 double type 															//memory requirement, hence 2*size extra memeory allocated

	if (!intSendPack || !intAllRecvPack || !doubleSendPack || !doubleAllRecvPack) 
	{
		printf("ERROR in prioritize_move: process:%d does not have sufficient memory, process exiting...\n", rank);
		exit(1);
	}

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;
	int movingIndex = 0;

	int elementCount = 0;

	srand(time(0));
	

	for (int j = 0; j < myPortionVertexSize; j++) 
	{
		int target = rand() % myPortionVertexSize;
		// swap numbers between i and target.
		int tmp = myVertices[j];
		myVertices[j] = myVertices[target];
		myVertices[target] = tmp;
	}

	//mark 1 end

	//__itt_resume();

	omp_set_num_threads(numTh);

#pragma omp parallel for
	for (int i = 0; i < myPortionVertexSize; i++) 
	{
		int vertexIndex = myVertices[i];
		Node& nd = nodes[vertexIndex];
		int oldMod = nd.ModIdx();

		int tid = omp_get_thread_num ();
			
		int nModLinks = 0;

		flowmap outFlowToMod; // <modId, flow> for outFlow...
		flowmap inFlowFromMod; // <modId, flow> for inFlow...
		
		flowmap::iterator iter;	
			
		auto tik = std::chrono::high_resolution_clock::now();

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) 
		{
			int newMod = nodes[linkIt->first].ModIdx();
			
			iter = outFlowToMod.find(newMod);

			if (iter != outFlowToMod.end()) 
			{
				iter->second += beta * linkIt->second;
			} 
			else 
			{
				outFlowToMod[newMod] = beta * linkIt->second; // initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;				
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) 
		{
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) 
			{
				inFlowFromMod[newMod] += beta * linkIt->second;
			} 
			else 
			{
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size()) {
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::prioritize_move()."
					<< endl;
		}

		double ndSize = nd.Size();
		double ndTPWeight = nd.TeleportWeight();
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;

		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize + beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize) + beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		int newMod;

		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) 
		{
			newMod = it->first;

			if (newMod == oldMod) 
			{
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} 
			else 
			{
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr + beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;

		if (outFlowToMod.count(oldMod) > 0) 
		{
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node

		bestResult.newModule = oldMod;

		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) 
		{
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) 
			{
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);

				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);

				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) 
				{

					bestResult.newModule = currentResult.newModule;

					bestResult.diffCodeLen = currentResult.diffCodeLen;

					bestResult.sumPr1 = currentResult.sumPr1;

					bestResult.sumPr2 = currentResult.sumPr2;

					bestResult.exitPr1 = currentResult.exitPr1;

					bestResult.exitPr2 = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;
				}
			}
		}
		auto ding = std::chrono::high_resolution_clock::now();
		if (bestResult.diffCodeLen < vThresh) 
		{
			int newMod = bestResult.newModule;

			if (oldMod < newMod) 
			{
				#pragma omp critical		
				{
					if (modules[newMod].numMembers == 0) 
					{
						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					modules[newMod].numMembers++;

					double exitPrChangeNewMod = bestResult.exitPr2 - modules[newMod].exitPr;

					modules[newMod].exitPr = bestResult.exitPr2;

					modules[newMod].sumPr = bestResult.sumPr2;

					modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;

					modules[newMod].sumTPWeight += ndTPWeight;

					if (nd.IsDangling()) 
					{
						modules[newMod].sumDangling += ndSize;
						modules[oldMod].sumDangling -= ndSize;
					}

					modules[oldMod].numMembers--;

					double exitPrChangeOldMod = modules[oldMod].exitPr - bestResult.exitPr1;

					modules[oldMod].exitPr = bestResult.exitPr1;

					modules[oldMod].sumPr = bestResult.sumPr1;

					modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;

					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) 
					{
						nEmptyMod++;
						nModule--;
					}

					intSendPack[movingIndex] = vertexIndex;

					intSendPack[movingIndex + 1] = bestResult.newModule;

					doubleSendPack[movingIndex] = exitPrChangeOldMod;

					doubleSendPack[movingIndex + 1] = exitPrChangeNewMod;

					sumAllExitPr = bestResult.newSumExitPr;

					codeLengthReduction += bestResult.diffCodeLen;

					numMoved++;
					movingIndex += 2;
				}
				//  This can be done in parallel without any locking, since the written value is always same.

				for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) 
				{

					isActives[linkIt->first] = 1;	// set as an active nodes.
				}
				for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) 
				{

					isActives[linkIt->first] = 1;	// set as an active nodes.
				}
			} 
		}
		auto dong = std::chrono::high_resolution_clock::now();

		if(tid == 13)
		bestModuleSelectionTime += std::chrono::duration_cast<std::chrono::nanoseconds>(dong - ding).count();

		#pragma omp critical
		{
			elementCount++;
		}//END critical

		auto tok = std::chrono::high_resolution_clock::now();
	
		if(tid == 13)
		prioritizeMoveNodeAccessTimeChunk += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

	}

	//mark 5 begin

	//__itt_pause();	

	if (myPortionVertexSize != elementCount) {
		printf(
				"something is not right, for rank:%d, myPortionVertexSize:%d does not match elementCount:%d\n",
				rank, myPortionVertexSize, elementCount);
	}

	
	MPI_Allgather(&numMoved, 1, MPI_INT, sendableElementCounts, 1, MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < totalProcessor; i++) 
	{
		intCounts[i] = 2 * sendableElementCounts[i];
		doubleCounts[i] = intCounts[i] + 2; //for storing sumAllExitPr and codelength in every process
	}

	intDispls[0] = 0;
	doubleDispls[0] = 0;

	for (int i = 1; i < totalProcessor; i++) 
	{
		intDispls[i] = intDispls[i - 1] + intCounts[i - 1];
		doubleDispls[i] = doubleDispls[i - 1] + doubleCounts[i - 1];
	}

	doubleSendPack[movingIndex] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[movingIndex + 1] = codeLengthReduction;

	auto ting = std::chrono::high_resolution_clock::now();

	MPI_Allgatherv(intSendPack, intCounts[rank], MPI_INT, intAllRecvPack,
			intCounts, intDispls, MPI_INT, MPI_COMM_WORLD);

	MPI_Allgatherv(doubleSendPack, doubleCounts[rank], MPI_DOUBLE,
			doubleAllRecvPack, doubleCounts, doubleDispls, MPI_DOUBLE,
			MPI_COMM_WORLD);

	auto tong = std::chrono::high_resolution_clock::now();

	prioritizeNodeCommunicationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tong - ting).count();

	auto ping = std::chrono::high_resolution_clock::now();	

	for (int p = 0; p < totalProcessor; p++) 
	{

		if (p != rank) 
		{
			int startIntIndex = intDispls[p];
			int endIntIndex = intDispls[p] + intCounts[p];

			//int totalElementSentFromSender = intAllRecvPack[endIntIndex];

			int totalElementSentFromSender = sendableElementCounts[p];

			int startDblIndex = doubleDispls[p];
			int endDblIndex = doubleDispls[p] + doubleCounts[p];

			sumAllExitPr += doubleAllRecvPack[endDblIndex - 2];
			codeLength += doubleAllRecvPack[endDblIndex - 1];

			//int highest_index = totalElementSentFromSender + startIntIndex;

			for (int z = startIntIndex, d = startDblIndex; z < endIntIndex; z += 2, d += 2) 
			{

				int nodeIndex = intAllRecvPack[z];

				Node& nod = nodes[nodeIndex];
				int oldModule = nod.modIdx;

				int newModule = intAllRecvPack[z + 1];

				nod.setModIdx(newModule);

				if (oldModule < newModule) 
				{
					numMoved++;

					if (modules[newModule].numMembers == 0) 
					{
						nEmptyMod--;
						nModule++;
					}

					modules[newModule].numMembers++;
					modules[newModule].exitPr += doubleAllRecvPack[d + 1]; //now how to solve that alternating problems??

					modules[newModule].sumPr = modules[newModule].sumPr
							+ nod.size;

					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;

					modules[newModule].sumTPWeight += nod.TeleportWeight();

					if (nod.IsDangling()) {
						modules[newModule].sumDangling += nod.Size();
						modules[oldModule].sumDangling -= nod.Size();
					}

					modules[oldModule].numMembers--;

					modules[oldModule].exitPr -= doubleAllRecvPack[d];

					if (modules[oldModule].exitPr < 0.0) 
					{
						modules[oldModule].exitPr = 0.0;
					}

					modules[oldModule].sumPr = modules[oldModule].sumPr
							- nod.size;

					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;

					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					if (modules[oldModule].numMembers == 0) 
					{
						nEmptyMod++;
						nModule--;
					}

					// update activeNodes and isActives vectors.
					// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
					for (link_iterator linkIt = nod.outLinks.begin(); linkIt != nod.outLinks.end(); linkIt++) 
					{
						isActives[linkIt->first] = 1;// set as an active nodes.
					}
					for (link_iterator linkIt = nod.inLinks.begin(); linkIt != nod.inLinks.end(); linkIt++) 
					{
						isActives[linkIt->first] = 1;// set as an active nodes.
					}
				}
			}
		}
	}

	vector<int>().swap(activeNodes);
	map<int, int>().swap(activeVertices);

	for (int i = 0; i < isActives.size(); i++) 
	{
		if (isActives[i] == 1) 
		{
			activeNodes.push_back(i);
			activeVertices.insert(make_pair(i, 1));
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	auto pong = std::chrono::high_resolution_clock::now();

	prioritizeNodeSyncTime += std::chrono::duration_cast<std::chrono::nanoseconds>(pong - ping).count();


	delete[] counts;
	delete[] intCounts;
	delete[] intDispls;
	delete[] doubleCounts;
	delete[] doubleDispls;
	delete[] intSendPack;
	delete[] intAllRecvPack;
	delete[] doubleSendPack;
	delete[] doubleAllRecvPack;


	auto end_time = std::chrono::high_resolution_clock::now();

	prioritizeMoveTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
	
	return numMoved;

}


/*
 *	+ This is a parallel implementation of Network::move() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */
int Network::parallelMove(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

// Generate random sequential order of nodes.
	vector<int> randomOrder(nNode);
	for (int i = 0; i < nNode; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nNode; i++) {
		int target = R->randInt(nNode - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);
	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nNode; i++) {
		Node& nd = nodes[randomOrder[i]];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		// Below values are related current module information, but it could be changed in the middle of this decision process.
		// However, we used this snapshot for finding next module of current node.
		// These will be a guideline of decision process and the correct values will be calculated at the end of this iteration.
		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		if (bestResult.diffCodeLen < 0.0) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					modules[newMod].numMembers++;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (nd.IsDangling()) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers--;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved++;	// This is inside of critical section, so no race condition happens.
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

		}

	}	// END parallel for

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

/*
 *	+ This is a parallel implementation of Network::prioritize_move() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */
int Network::prioritize_parallelMove(int numTh, double& tSequential,
		double vThresh) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nActive = activeNodes.size();
	int nNextActive = 0;

// Generate random sequential order of nodes.
	vector<int> randomOrder(nActive);
	for (int i = 0; i < nActive; i++)
		randomOrder[i] = activeNodes[i];

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

// Now randomOrder vector already had all activeNodes info,
// so we can reset activeNodes vector for adding for new active nodes.
	vector<int>().swap(activeNodes);

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);
	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nActive; i++) {
		Node& nd = nodes[randomOrder[i]];// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		int nModLinks = 0;// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = nodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout
					<< "ALERT: nModLinks != outFlowToMod.size() in Network::move()."
					<< endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		// Below values are related current module information, but it could be changed in the middle of this decision process.
		// However, we used this snapshot for finding next module of current node.
		// These will be a guideline of decision process and the correct values will be calculated at the end of this iteration.
		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//////////////////// THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].members.size() > 1 && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}
		}

		// store the best possilbe movement information if necessary.
		// In this version, we are trying to move as soon as possible, i.e. right after decision making...
		//if (bestResult.diffCodeLen < 0.0) {
		if (bestResult.diffCodeLen < vThresh) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

				// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = nodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = nodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					modules[newMod].numMembers++;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (nd.IsDangling()) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers--;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved++;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

			// update activeNodes and isActives vectors.
			// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
			// This can be done in parallel without any locking, since the written value is always same.
			for (link_iterator linkIt = nd.outLinks.begin();
					linkIt != nd.outLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.

			for (link_iterator linkIt = nd.inLinks.begin();
					linkIt != nd.inLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;		// set as an active nodes.
		}

	}	// END parallel for

	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	if (nodes.size() != nModule + nEmptyMod) {
		cout << "Something wrong!! nodes.size() != nModule + nEmptyMod."
				<< endl;
	}

	return numMoved;
}

/*
 *	same as Network::move() except that the moving unit is SuperNode.
 */

int Network::moveSuperNodes(int iteration) {

	int nSuperNodes = superNodes.size();

	int nNextActive = 0;
	int SPNodeCountforSending = 0;
	int elementCount = 0;

	bool incrementElementCount = false;
	set<int> emptyModuleSet;

	int size;
	int rank;

	struct timeval startMoveSuperNodes, endMoveSuperNodes;

	gettimeofday(&startMoveSuperNodes, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// Generate random sequential order of nodes.
	int* randomGlobalArray = new int[nSuperNodes]();

	int start, end;
	findAssignedPart(&start, &end, nSuperNodes, size, rank);

	int numberOfSPNode = end - start;//number of elements assigned to this processor

	for (int i = start; i < end; i++) {
		randomGlobalArray[i] = i;
	}

	for (int i = start; i < end; i++) {
		int pretarget = rand() % numberOfSPNode;
		int target = start + pretarget;
		// swap numbers between i and target.
		int tmp = randomGlobalArray[i];
		randomGlobalArray[i] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}

	double totalElemenets = nSuperNodes;

	int stripSize = ceil(totalElemenets / size) + 1;

	const unsigned int intPackSize = 2 * stripSize + 3; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 4 * stripSize + 3;

	int* intSendPack = new (nothrow) int[intPackSize](); // multiplier 3 for 3 different elements of index, newModules, oldModules and +1 for elementCount
	int* intReceivePack = new (nothrow) int[intPackSize]();
	int* intAllRecvPack = new (nothrow) int[intPackSize * size]();

	double* doubleSendPack = new (nothrow) double[doublePackSize](); // multiplier 6 for 6 different elements of diffCodeLen, sumPr1, sumPr2, exitPr1, exitPr2, newSumExitPr
	double* doubleReceivePack = new (nothrow) double[doublePackSize]();

	if (!intSendPack || !intReceivePack || !doubleSendPack
			|| !doubleReceivePack) {
		printf(
				"ERROR in moveSuperNodes: process:%d does not have sufficient memory, process exiting...\n",
				rank);
		exit(1);
	}

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;

	SPNodeCountforSending = 0;

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = start; i < end; i++) {

		int superNodeIndex = randomGlobalArray[i];
		SuperNode& nd = superNodes[superNodeIndex]; // look at i_th Node of the random sequential order.

		bool moduleUpdate = false;

		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0; // The number of links to/from between the current node and other modules.

		flowmap outFlowToMod; // <modID, flow> for outFlow...
		flowmap inFlowFromMod; // <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second; // initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size(); // p_nd.
		double ndTPWeight = nd.TeleportWeight(); // tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		//vector<MoveSummary> moveResults(nModLinks);
		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod; // I have to do that for the sole purpose of sending the information of moving to a new module when sending through MPI although may be there is no module change at all

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize; // This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {

					//save new module id
					bestResult.newModule = currentResult.newModule;

					bestResult.diffCodeLen = currentResult.diffCodeLen;

					bestResult.sumPr1 = currentResult.sumPr1;

					bestResult.sumPr2 = currentResult.sumPr2;

					bestResult.exitPr1 = currentResult.exitPr1;

					bestResult.exitPr2 = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			}
		}

		// Make best possible move for the current node nd.
		if (bestResult.diffCodeLen < 0.0) {
			// update related to newMod...
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			if (oldMod < newMod) {

				if (modules[newMod].numMembers == 0) {
					nEmptyMod--;
					nModule++;
				}

				nd.setModIdx(newMod);

				for (int j = 0; j < spMembers; j++) {
					nd.members[j]->setModIdx(newMod);
				}

				modules[newMod].numMembers += spMembers;
				modules[newMod].exitPr = bestResult.exitPr2;
				modules[newMod].sumPr = bestResult.sumPr2;
				modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
				modules[newMod].sumTPWeight += ndTPWeight;

				double spDanglingSize = nd.DanglingSize();
				if (spDanglingSize > 0.0) {
					modules[newMod].sumDangling += spDanglingSize;
					modules[oldMod].sumDangling -= spDanglingSize;
				}

				// update related to the oldMod...
				modules[oldMod].numMembers -= spMembers;
				modules[oldMod].exitPr = bestResult.exitPr1;
				modules[oldMod].sumPr = bestResult.sumPr1;
				modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
				modules[oldMod].sumTPWeight -= ndTPWeight;

				if (modules[oldMod].numMembers == 0) {
					nEmptyMod++;
					nModule--;
				}

				//the following code is for sending module update across processors
				intSendPack[SPNodeCountforSending] = superNodeIndex;
				intSendPack[1 * stripSize + SPNodeCountforSending] =
						bestResult.newModule;

				doubleSendPack[SPNodeCountforSending] = bestResult.sumPr1;

				doubleSendPack[1 * stripSize + SPNodeCountforSending] =
						bestResult.sumPr2;

				doubleSendPack[2 * stripSize + SPNodeCountforSending] =
						bestResult.exitPr1;

				doubleSendPack[3 * stripSize + SPNodeCountforSending] =
						bestResult.exitPr2;

				sumAllExitPr = bestResult.newSumExitPr;

				SPNodeCountforSending++;

				codeLengthReduction += bestResult.diffCodeLen;
				numMoved += spMembers; // although we moved a superNode, a superNode is actually a set of nodes.
			}
		}

		elementCount++;

	}

	if (numberOfSPNode != elementCount) {
		printf(
				"something is not right in prioritize_movespnodes, for rank:%d, numberOfSPNode:%d does not match elementCount:%d\n",
				rank, numberOfSPNode, elementCount);
	}

	intSendPack[intPackSize - 1] = SPNodeCountforSending;
	doubleSendPack[doublePackSize - 1] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[doublePackSize - 2] = codeLengthReduction;

	/*	MPI_Allgather(intSendPack, intPackSize, MPI_INT, intAllRecvPack,
	 intPackSize,
	 MPI_INT, MPI_COMM_WORLD);*/

	for (int processId = 0; processId < size; processId++) {
		if (processId != rank) {
			MPI_Sendrecv(intSendPack, intPackSize, MPI_INT, processId,
					iteration, intReceivePack, intPackSize, MPI_INT, processId,
					iteration,
					MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			MPI_Sendrecv(doubleSendPack, doublePackSize, MPI_DOUBLE, processId,
					iteration, doubleReceivePack, doublePackSize,
					MPI_DOUBLE, processId, iteration, MPI_COMM_WORLD,
					MPI_STATUSES_IGNORE);

			int totalSPNodefromSender = intReceivePack[intPackSize - 1];

			sumAllExitPr += doubleReceivePack[doublePackSize - 1];
			codeLength += doubleReceivePack[doublePackSize - 2];

			for (int z = 0; z < totalSPNodefromSender; z++) {
				int nodeIndex = intReceivePack[z];
				SuperNode& nod = superNodes[nodeIndex]; // I didn't need to send the indices from the processes because every process has already this information from counts and displs

				int spMembers = nod.members.size();

				int oldModule = nod.modIdx;

				int newModule = intReceivePack[1 * stripSize + z];

				if (oldModule != newModule) {

					numMoved += spMembers;

					if (modules[newModule].numMembers == 0) {
						nEmptyMod--;
						nModule++;
					}

					nod.setModIdx(newModule);

					for (int j = 0; j < spMembers; j++) {
						nod.members[j]->setModIdx(newModule);
					}

					modules[newModule].numMembers += spMembers;
					modules[newModule].exitPr = doubleReceivePack[3 * stripSize
							+ z];
					modules[newModule].sumPr = doubleReceivePack[1 * stripSize
							+ z];
					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;
					modules[newModule].sumTPWeight += nod.TeleportWeight();

					double spDanglingSize = nod.DanglingSize();

					if (spDanglingSize > 0.0) {
						modules[newModule].sumDangling += spDanglingSize;
						modules[oldModule].sumDangling -= spDanglingSize;
					}

					modules[oldModule].numMembers -= spMembers;
					modules[oldModule].exitPr = doubleReceivePack[2 * stripSize
							+ z];
					modules[oldModule].sumPr = doubleReceivePack[z];
					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;
					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					if (modules[oldModule].numMembers == 0) {
						nEmptyMod++;
						nModule--;
					}
				}
			}
		}
	}

	delete[] randomGlobalArray;
	delete[] intSendPack;
	delete[] intReceivePack;
	delete[] doubleSendPack;
	delete[] doubleReceivePack;

	gettimeofday(&endMoveSuperNodes, NULL);

	return numMoved;
}

void Network::showOutput(int iteration, int prioritizeSPMove, bool inWhile) {
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < nNode; i++) {
		printf(
				"before distributed output for rank:%d, iteration:%d, nodes:%d, module Id:%d\n",
				rank, iteration, i, nodes[i].modIdx);
	}
	//printf("before distributed output for rank:%d, number of superNodes:%d\n",
			//rank, superNodes.size());
}

/*
 * This function implements a prioritized version of Network::moveSuperNodes() above.
 */

int Network::prioritize_moveSPnodes(double vThresh, int numTh, int iteration, bool inWhile) 
{

	int SPNodeCountforSending = 0;
	int elementCount = 0;

	srand(time(0));

	int size;
	int rank;

	auto start_time = std::chrono::high_resolution_clock::now();

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	iteration_prioritizeSpMove++;

	int nActive = activeNodes.size();

	double totalElements = nActive;

	int stripSize = ceil(totalElements / size) + 1;

	int* randomGlobalArray = new int[nActive]();

	int start, end;

	findAssignedPart(&start, &end, nActive, size, rank);

	int numberOfSPNode = end - start; //number of elements assigned to this processor

	for (int i = start; i < end; i++) 
	{
		randomGlobalArray[i] = activeNodes[i];
	}	

	for (int j = start; j < end; j++) 
	{
		int preTarget = rand() % numberOfSPNode;
		int target = start + preTarget;
		// swap numbers between i and target.
		int tmp = randomGlobalArray[j];
		randomGlobalArray[j] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}

	const int totalProcessor = size;

	int *sendableElementCounts = new int[totalProcessor]();

	int *intCounts = new int[totalProcessor]();
	int *intDispls = new int[totalProcessor]();

	int *doubleCounts = new int[totalProcessor]();
	int *doubleDispls = new int[totalProcessor]();

// Generate random sequential order of active nodes.
	const unsigned int intPackSize = 2 * stripSize; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 2 * stripSize;

	int* intSendPack = new (nothrow) int[intPackSize]();
	//int* intReceivePack = new (nothrow) int[intPackSize]();
	int* intAllRecvPack = new (nothrow) int[intPackSize * size]();

	double* doubleSendPack = new (nothrow) double[doublePackSize]();
	//double* doubleReceivePack = new (nothrow) double[doublePackSize]();
	double* doubleAllRecvPack = new (nothrow) double[doublePackSize * size]();

	if (!intSendPack || !intAllRecvPack || !doubleSendPack || !doubleAllRecvPack) 
	{
		printf("ERROR in prioritize_moveSPnodes: process:%d does not have sufficient memory, process exiting...\n", rank);
		exit(1);
	}

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;
	int movingIndex = 0;

	SPNodeCountforSending = 0;
	
	//__itt_resume();
	
	omp_set_num_threads(numTh);
	
// Move each node to one of its neighbor modules in random sequential order.
#pragma omp parallel for
	for (int i = start; i < end; i++) 
	{

		int superNodeIndex = randomGlobalArray[i];
		SuperNode& nd = superNodes[superNodeIndex]; // look at i_th Node of the random sequential order.

		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0; // The number of links to/from between the current node and other modules.

		flowmap outFlowToMod; // <modID, flow> for outFlow...
		flowmap inFlowFromMod; // <modID, flow> for inFlow...

		// count other modules that are connected to the current node.

		int tid = omp_get_thread_num ();

		auto tik = std::chrono::high_resolution_clock::now();

		for (link_iterator linkIt = nd.outLinks.begin(); linkIt != nd.outLinks.end(); linkIt++) 
		{
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) 
			{
				outFlowToMod[newMod] += beta * linkIt->second;
			} 
			else 
			{
				outFlowToMod[newMod] = beta * linkIt->second; // initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin(); linkIt != nd.inLinks.end(); linkIt++) 
		{
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) 
			{
				inFlowFromMod[newMod] += beta * linkIt->second;
			} 
			else 
			{
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size(); // p_nd.
		double ndTPWeight = nd.TeleportWeight(); // tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.

		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) 
		{
			int newMod = it->first;
			if (newMod == oldMod) 
			{
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} 
			else 
			{
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize) * modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr + beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) 
		{
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod + inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod; // I have to do that for the sole purpose of sending the information of moving to a new module when sending through MPI although may be there is no module change at all

		for (flowmap::iterator it = outFlowToMod.begin(); it != outFlowToMod.end(); it++) 
		{
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) 
			{
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize; // This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr() - outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1 + currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1) + pLogP(currentResult.exitPr2) - pLogP(oldExitPr1) - pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(currentResult.exitPr1 + currentResult.sumPr1) + pLogP(currentResult.exitPr2 + currentResult.sumPr2) - pLogP(oldExitPr1 + oldSumPr1) - pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) 
				{
					//save new module id
					bestResult.newModule = currentResult.newModule;

					bestResult.diffCodeLen = currentResult.diffCodeLen;

					bestResult.sumPr1 = currentResult.sumPr1;

					bestResult.sumPr2 = currentResult.sumPr2;

					bestResult.exitPr1 = currentResult.exitPr1;

					bestResult.exitPr2 = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			}
		}

		// Make best possible move for the current node nd.
		//if (bestResult.diffCodeLen < 0.0) {

		auto ding = std::chrono::high_resolution_clock::now();
		if (bestResult.diffCodeLen < vThresh) 
		{
			// update related to newMod...

			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			if (oldMod < newMod) 
			{
				
				#pragma omp critical		
				{
					if (modules[newMod].numMembers == 0) 
					{
						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);

					for (int j = 0; j < spMembers; j++) 
					{
						nd.members[j]->setModIdx(newMod);
					}

					modules[newMod].numMembers += spMembers;

					double changeNewModule = bestResult.exitPr2 - modules[newMod].exitPr;

					modules[newMod].exitPr = bestResult.exitPr2;

					modules[newMod].sumPr = bestResult.sumPr2;
					modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					double spDanglingSize = nd.DanglingSize();
					if (spDanglingSize > 0.0) 
					{
						modules[newMod].sumDangling += spDanglingSize;
						modules[oldMod].sumDangling -= spDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers -= spMembers;

					double changeOldModule = modules[oldMod].exitPr - bestResult.exitPr1;

					modules[oldMod].exitPr = bestResult.exitPr1;

					modules[oldMod].sumPr = bestResult.sumPr1;
					modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) 
					{
						nEmptyMod++;
						nModule--;
					}

					//the following code is for sending module update across processors
					intSendPack[movingIndex] = superNodeIndex;

					intSendPack[1 + movingIndex] = bestResult.newModule;

					doubleSendPack[movingIndex] = changeOldModule;

					doubleSendPack[1 + movingIndex] = changeNewModule;

					sumAllExitPr = bestResult.newSumExitPr;

					codeLengthReduction += bestResult.diffCodeLen;

					SPNodeCountforSending++;

					movingIndex += 2;

					numMoved += spMembers;
				}
				// update activeNodes and isActives vectors.
				// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
				for (link_iterator linkIt = nd.outLinks.begin();
						linkIt != nd.outLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.

				for (link_iterator linkIt = nd.inLinks.begin();
						linkIt != nd.inLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.*/
			}
		}
		
		auto dong = std::chrono::high_resolution_clock::now();

		if(tid == 13)
		bestSpModuleSelectionTime += std::chrono::duration_cast<std::chrono::nanoseconds>(dong - ding).count();
		
		#pragma omp critical
		{
			elementCount++;
		}

		auto tok = std::chrono::high_resolution_clock::now();
		if(tid == 13)
		prioritizeSpNodeModuleAccessTimeChunk += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

		threadtimes[tid] += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();
	}

	//__itt_pause();

	if (numberOfSPNode != elementCount) {
		if (iteration == 0) {
			printf(
					"something is not right in prioritize_movespnodes, for rank:%d, numberOfSPNode:%d does not match elementCount:%d\n",
					rank, numberOfSPNode, elementCount);
		}
	}


//gettimeofday(&startPrio_moveSPnodes_MPISendRecv, NULL);
	MPI_Allgather(&SPNodeCountforSending, 1, MPI_INT, sendableElementCounts, 1,
	MPI_INT, MPI_COMM_WORLD);

	for (int i = 0; i < totalProcessor; i++) 
	{
		intCounts[i] = 2 * sendableElementCounts[i];
		doubleCounts[i] = intCounts[i] + 2;
	}

	intDispls[0] = 0;
	doubleDispls[0] = 0;

	for (int i = 1; i < totalProcessor; i++) 
	{
		intDispls[i] = intDispls[i - 1] + intCounts[i - 1];
		doubleDispls[i] = doubleDispls[i - 1] + doubleCounts[i - 1];
	}

	doubleSendPack[movingIndex] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[movingIndex + 1] = codeLengthReduction;

	auto ting = std::chrono::high_resolution_clock::now();

	MPI_Allgatherv(intSendPack, intCounts[rank], MPI_INT, intAllRecvPack,
			intCounts, intDispls,
			MPI_INT, MPI_COMM_WORLD);

	MPI_Allgatherv(doubleSendPack, doubleCounts[rank], MPI_DOUBLE,
			doubleAllRecvPack, doubleCounts, doubleDispls, MPI_DOUBLE,
			MPI_COMM_WORLD);

	auto tong = std::chrono::high_resolution_clock::now();

	prioritizeSpNodeCommunicationTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tong - ting).count();	

	auto ping = std::chrono::high_resolution_clock::now();

	for (int p = 0; p < size; p++) 
	{
		if (p != rank) 
		{
			int startIntIndex = intDispls[p];
			int endIntIndex = intDispls[p] + intCounts[p];
			int totalElementSentFromSender = sendableElementCounts[p];
			int startDblIndex = doubleDispls[p];
			int endDblIndex = doubleDispls[p] + doubleCounts[p];

			codeLength += doubleAllRecvPack[endDblIndex - 1];
			sumAllExitPr += doubleAllRecvPack[endDblIndex - 2];

			for (int z = startIntIndex, d = startDblIndex; z < endIntIndex; z +=
					2, d += 2) {

				int nodeIndex = intAllRecvPack[z];

				SuperNode& nod = superNodes[nodeIndex];
				int spMembers = nod.members.size();
				int oldModule = nod.modIdx;
				int newModule = intAllRecvPack[1 + z];

				if (oldModule < newModule) {

					numMoved += spMembers;

					if (modules[newModule].numMembers == 0) {
						nEmptyMod--;
						nModule++;
					}

					nod.setModIdx(newModule);

					for (int j = 0; j < spMembers; j++) {
						nod.members[j]->setModIdx(newModule);
					}

					modules[newModule].numMembers += spMembers;
					modules[newModule].exitPr += doubleAllRecvPack[1 + d];
					/*					modules[newModule].sumPr = doubleAllRecvPack[1 * stripSize
					 + d];*/

					modules[newModule].sumPr = modules[newModule].sumPr
							+ nod.size;

					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;

					modules[newModule].sumTPWeight += nod.TeleportWeight();

					double spDanglingSize = nod.DanglingSize();

					if (spDanglingSize > 0.0) 
					{
						modules[newModule].sumDangling += spDanglingSize;
						modules[oldModule].sumDangling -= spDanglingSize;
					}

					modules[oldModule].numMembers -= spMembers;

					modules[oldModule].exitPr -= doubleAllRecvPack[d];

					modules[oldModule].sumPr = modules[oldModule].sumPr
							- nod.size;

					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;
					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					if (modules[oldModule].numMembers == 0) 
					{
						nEmptyMod++;
						nModule--;
					}

					// update activeNodes and isActives vectors.
					// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
					for (link_iterator linkIt = nod.outLinks.begin();
							linkIt != nod.outLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.

					for (link_iterator linkIt = nod.inLinks.begin();
							linkIt != nod.inLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.
				}
			}
		}
	}

	vector<int>().swap(activeNodes);
	for (int i = 0; i < isActives.size(); i++) 
	{
		if (isActives[i] == 1) 
		{
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	auto pong = std::chrono::high_resolution_clock::now();

	prioritizeSpNodeSyncTime += std::chrono::duration_cast<std::chrono::nanoseconds>(pong - ping).count();

	delete[] intCounts;
	delete[] intDispls;
	delete[] doubleCounts;
	delete[] doubleDispls;
	delete[] randomGlobalArray;
	delete[] intSendPack;
	delete[] intAllRecvPack;
	delete[] doubleSendPack;
	delete[] doubleAllRecvPack;

	auto end_time = std::chrono::high_resolution_clock::now();

	prioritizeSpNodeMoveTime += std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
	
	return numMoved;
}

/*
 * This function implements a prioritized version of Network::moveSuperNodes() above.
 */

int Network::old_prioritize_moveSPnodes(double vThresh, int tag, int iteration,
		bool inWhile, double& total_time_prioritize_Spmove,
		int& total_iterations_priorSPNodes, double& total_time_MPISendRecvSP) {

	int SPNodeCountforSending = 0;
	int elementCount = 0;
	int mpi_tag = 5555 + iteration;

	srand(time(0));

	set<int> emptyModuleSet;

	int size;
	int rank;

	/*struct timeval startPrio_moveSPnodes, endPrio_moveSPnodes;
	 struct timeval startPrio_moveSPnodes_MPISendRecv,
	 endPrio_moveSPnodes_MPISendRecv;*/

	double startPrio_moveSPnodes, endPrio_moveSPnodes;
	double startPrio_moveSPnodes_MPISendRecv, endPrio_moveSPnodes_MPISendRecv;

//gettimeofday(&startPrio_moveSPnodes, NULL);

	startPrio_moveSPnodes = MPI_Wtime();

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int nActive = activeNodes.size();

	double totalElements = nActive;

	int stripSize = ceil(totalElements / size) + 1;

	int* randomGlobalArray = new int[nActive]();

	int start, end;

	findAssignedPart(&start, &end, nActive, size, rank);

	int numberOfSPNode = end - start; //number of elements assigned to this processor

	for (int i = start; i < end; i++) {
		randomGlobalArray[i] = activeNodes[i];
	}

	for (int j = start; j < end; j++) {
		int preTarget = rand() % numberOfSPNode;
		int target = start + preTarget;
		// swap numbers between i and target.
		int tmp = randomGlobalArray[j];
		randomGlobalArray[j] = randomGlobalArray[target];
		randomGlobalArray[target] = tmp;
	}

// Generate random sequential order of active nodes.
	const unsigned int intPackSize = 2 * stripSize + 3; // 3*nActive+1 index will save emptyModCount and 3*nActive+2 index will save nModuleCount
	const unsigned int doublePackSize = 2 * stripSize + 3;

	int* intSendPack = new (nothrow) int[intPackSize]();
	//int* intReceivePack = new (nothrow) int[intPackSize]();
	int* intAllRecvPack = new (nothrow) int[intPackSize * size]();

	double* doubleSendPack = new (nothrow) double[doublePackSize]();
	//double* doubleReceivePack = new (nothrow) double[doublePackSize]();
	double* doubleAllRecvPack = new (nothrow) double[doublePackSize * size]();

	if (!intSendPack || !intAllRecvPack || !doubleSendPack
			|| !doubleAllRecvPack) {
		printf(
				"ERROR in prioritize_moveSPnodes: process:%d does not have sufficient memory, process exiting...\n",
				rank);
		exit(1);
	}

	double codeLengthReduction = 0.0;
	double currentSumAllExitPr = sumAllExitPr;

	int numMoved = 0;

	SPNodeCountforSending = 0;

// Move each node to one of its neighbor modules in random sequential order.
	for (int i = start; i < end; i++) {
		int superNodeIndex = randomGlobalArray[i];
		SuperNode& nd = superNodes[superNodeIndex]; // look at i_th Node of the random sequential order.

		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0; // The number of links to/from between the current node and other modules.

		flowmap outFlowToMod; // <modID, flow> for outFlow...
		flowmap inFlowFromMod; // <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second; // initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size(); // p_nd.
		double ndTPWeight = nd.TeleportWeight(); // tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

		// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0; // This is the default value, if we can't find diffCodeLen < 0, then don't move the node.
		bestResult.newModule = oldMod; // I have to do that for the sole purpose of sending the information of moving to a new module when sending through MPI although may be there is no module change at all

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {
				// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

				// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize; // This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

				// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

				// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					//save new module id
					bestResult.newModule = currentResult.newModule;

					bestResult.diffCodeLen = currentResult.diffCodeLen;

					bestResult.sumPr1 = currentResult.sumPr1;

					bestResult.sumPr2 = currentResult.sumPr2;

					bestResult.exitPr1 = currentResult.exitPr1;

					bestResult.exitPr2 = currentResult.exitPr2;

					bestResult.newSumExitPr = currentResult.newSumExitPr;

				}
			}
		}

		// Make best possible move for the current node nd.
		//if (bestResult.diffCodeLen < 0.0) {

		if (bestResult.diffCodeLen < vThresh) {
			// update related to newMod...

			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

			if (oldMod < newMod) {

				if (modules[newMod].numMembers == 0) {
					nEmptyMod--;
					nModule++;
				}

				nd.setModIdx(newMod);

				for (int j = 0; j < spMembers; j++) {
					nd.members[j]->setModIdx(newMod);
				}

				modules[newMod].numMembers += spMembers;

				double changeNewModule = bestResult.exitPr2
						- modules[newMod].exitPr;

				modules[newMod].exitPr = bestResult.exitPr2;

				modules[newMod].sumPr = bestResult.sumPr2;
				modules[newMod].stayPr = bestResult.exitPr2 + bestResult.sumPr2;
				modules[newMod].sumTPWeight += ndTPWeight;

				double spDanglingSize = nd.DanglingSize();
				if (spDanglingSize > 0.0) {
					modules[newMod].sumDangling += spDanglingSize;
					modules[oldMod].sumDangling -= spDanglingSize;
				}

				// update related to the oldMod...
				modules[oldMod].numMembers -= spMembers;

				double changeOldModule = modules[oldMod].exitPr
						- bestResult.exitPr1;

				modules[oldMod].exitPr = bestResult.exitPr1;

				modules[oldMod].sumPr = bestResult.sumPr1;
				modules[oldMod].stayPr = bestResult.exitPr1 + bestResult.sumPr1;
				modules[oldMod].sumTPWeight -= ndTPWeight;

				if (modules[oldMod].numMembers == 0) {
					nEmptyMod++;
					nModule--;
				}

				//the following code is for sending module update across processors
				intSendPack[SPNodeCountforSending] = superNodeIndex;

				intSendPack[1 * stripSize + SPNodeCountforSending] =
						bestResult.newModule;

				doubleSendPack[SPNodeCountforSending] = changeOldModule;

				doubleSendPack[1 * stripSize + SPNodeCountforSending] =
						changeNewModule;

				sumAllExitPr = bestResult.newSumExitPr;

				SPNodeCountforSending++;

				codeLengthReduction += bestResult.diffCodeLen;

				numMoved += spMembers;

				// update activeNodes and isActives vectors.
				// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
				for (link_iterator linkIt = nd.outLinks.begin();
						linkIt != nd.outLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.

				for (link_iterator linkIt = nd.inLinks.begin();
						linkIt != nd.inLinks.end(); linkIt++)
					isActives[linkIt->first] = 1;	// set as an active nodes.*/
			}
		}

		elementCount++;
	}

	if (numberOfSPNode != elementCount) {
		if (iteration == 0) {
			printf(
					"something is not right in prioritize_movespnodes, for rank:%d, numberOfSPNode:%d does not match elementCount:%d\n",
					rank, numberOfSPNode, elementCount);
		}
	}

	intSendPack[intPackSize - 1] = SPNodeCountforSending;
	doubleSendPack[doublePackSize - 1] = sumAllExitPr - currentSumAllExitPr;
	codeLength += codeLengthReduction;
	doubleSendPack[doublePackSize - 2] = codeLengthReduction;

//gettimeofday(&startPrio_moveSPnodes_MPISendRecv, NULL);

	startPrio_moveSPnodes_MPISendRecv = MPI_Wtime();

	MPI_Allgather(intSendPack, intPackSize, MPI_INT, intAllRecvPack,
			intPackSize,
			MPI_INT, MPI_COMM_WORLD);

	MPI_Allgather(doubleSendPack, doublePackSize, MPI_DOUBLE, doubleAllRecvPack,
			doublePackSize, MPI_DOUBLE, MPI_COMM_WORLD);

	endPrio_moveSPnodes_MPISendRecv = MPI_Wtime();

//gettimeofday(&endPrio_moveSPnodes_MPISendRecv, NULL);

	/*total_time_MPISendRecvSP += elapsedTimeInSec(
	 startPrio_moveSPnodes_MPISendRecv, endPrio_moveSPnodes_MPISendRecv);
	 */

	total_time_MPISendRecvSP += endPrio_moveSPnodes_MPISendRecv
			- startPrio_moveSPnodes_MPISendRecv;

	for (int p = 0; p < size; p++) {
		if (p != rank) {
			int startIntIndex = p * intPackSize;
			int endIntIndex = ((p + 1) * intPackSize) - 1;
			int totalElementSentFromSender = intAllRecvPack[endIntIndex];
			int startDblIndex = p * doublePackSize;
			int endDblIndex = (p + 1) * doublePackSize;
			sumAllExitPr += doubleAllRecvPack[endDblIndex - 1];
			codeLength += doubleAllRecvPack[endDblIndex - 2];
			int highest_index = totalElementSentFromSender + startIntIndex;

			for (int z = startIntIndex, d = startDblIndex; z < highest_index;
					z++, d++) {

				int nodeIndex = intAllRecvPack[z];

				SuperNode& nod = superNodes[nodeIndex];
				int spMembers = nod.members.size();
				int oldModule = nod.modIdx;
				int newModule = intAllRecvPack[1 * stripSize + z];

				if (oldModule < newModule) {

					numMoved += spMembers;

					if (modules[newModule].numMembers == 0) {
						nEmptyMod--;
						nModule++;
					}

					nod.setModIdx(newModule);

					for (int j = 0; j < spMembers; j++) {
						nod.members[j]->setModIdx(newModule);
					}

					modules[newModule].numMembers += spMembers;
					modules[newModule].exitPr += doubleAllRecvPack[1 * stripSize
							+ d];
					/*					modules[newModule].sumPr = doubleAllRecvPack[1 * stripSize
					 + d];*/

					modules[newModule].sumPr = modules[newModule].sumPr
							+ nod.size;

					modules[newModule].stayPr = modules[newModule].exitPr
							+ modules[newModule].sumPr;

					modules[newModule].sumTPWeight += nod.TeleportWeight();

					double spDanglingSize = nod.DanglingSize();

					if (spDanglingSize > 0.0) {
						modules[newModule].sumDangling += spDanglingSize;
						modules[oldModule].sumDangling -= spDanglingSize;
					}

					modules[oldModule].numMembers -= spMembers;

					modules[oldModule].exitPr -= doubleAllRecvPack[d];

					/*					if (modules[oldModule].exitPr < 0.0) {
					 modules[oldModule].exitPr = 0.0;
					 }*/

					/*modules[oldModule].sumPr = doubleAllRecvPack[d];*/

					modules[oldModule].sumPr = modules[oldModule].sumPr
							- nod.size;

					modules[oldModule].stayPr = modules[oldModule].exitPr
							+ modules[oldModule].sumPr;
					modules[oldModule].sumTPWeight -= nod.TeleportWeight();

					if (modules[oldModule].numMembers == 0) {
						nEmptyMod++;
						nModule--;
					}

					// update activeNodes and isActives vectors.
					// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
					for (link_iterator linkIt = nod.outLinks.begin();
							linkIt != nod.outLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.

					for (link_iterator linkIt = nod.inLinks.begin();
							linkIt != nod.inLinks.end(); linkIt++)
						isActives[linkIt->first] = 1;// set as an active nodes.
				}
			}
		}
	}
	/*
	 for (int processId = 0; processId < size; processId++) {
	 if (processId != rank) {

	 gettimeofday(&startPrio_moveSPnodes_MPISendRecv, NULL);

	 MPI_Sendrecv(intSendPack, intPackSize, MPI_INT, processId,
	 mpi_tag, intReceivePack, intPackSize,
	 MPI_INT, processId, mpi_tag,
	 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

	 MPI_Sendrecv(doubleSendPack, doublePackSize, MPI_DOUBLE,
	 processId, mpi_tag, doubleReceivePack, doublePackSize,
	 MPI_DOUBLE, processId, mpi_tag, MPI_COMM_WORLD,
	 MPI_STATUSES_IGNORE);

	 gettimeofday(&endPrio_moveSPnodes_MPISendRecv, NULL);

	 total_time_MPISendRecvSP += elapsedTimeInSec(
	 startPrio_moveSPnodes_MPISendRecv,
	 endPrio_moveSPnodes_MPISendRecv);

	 int totalSPNodefromSender = intReceivePack[intPackSize - 1];

	 sumAllExitPr += doubleReceivePack[doublePackSize - 1];
	 codeLength += doubleReceivePack[doublePackSize - 2];

	 for (int z = 0; z < totalSPNodefromSender; z++) {

	 int nodeIndex = intReceivePack[z];
	 SuperNode& nod = superNodes[nodeIndex];	// I didn't need to send the indices from the processes because every process has already this information from counts and displs

	 int spMembers = nod.members.size();

	 int oldModule = nod.modIdx;

	 int newModule = intReceivePack[1 * stripSize + z];

	 if (oldModule != newModule) {

	 numMoved += spMembers;

	 if (modules[newModule].numMembers == 0) {
	 nEmptyMod--;
	 nModule++;
	 }

	 nod.setModIdx(newModule);

	 for (int j = 0; j < spMembers; j++) {
	 nod.members[j]->setModIdx(newModule);
	 }

	 modules[newModule].numMembers += spMembers;

	 modules[newModule].exitPr = doubleReceivePack[3
	 * stripSize + z];
	 modules[newModule].sumPr = doubleReceivePack[1
	 * stripSize + z];
	 modules[newModule].stayPr = modules[newModule].exitPr
	 + modules[newModule].sumPr;
	 modules[newModule].sumTPWeight += nod.TeleportWeight();

	 double spDanglingSize = nod.DanglingSize();

	 if (spDanglingSize > 0.0) {
	 modules[newModule].sumDangling += spDanglingSize;
	 modules[oldModule].sumDangling -= spDanglingSize;
	 }

	 modules[oldModule].numMembers -= spMembers;

	 modules[oldModule].exitPr = doubleReceivePack[2
	 * stripSize + z];

	 modules[oldModule].sumPr = doubleReceivePack[z];

	 modules[oldModule].stayPr = modules[oldModule].exitPr
	 + modules[oldModule].sumPr;
	 modules[oldModule].sumTPWeight -= nod.TeleportWeight();

	 if (modules[oldModule].numMembers == 0) {
	 nEmptyMod++;
	 nModule--;
	 }

	 // update activeNodes and isActives vectors.
	 // We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
	 for (link_iterator linkIt = nod.outLinks.begin();
	 linkIt != nod.outLinks.end(); linkIt++)
	 isActives[linkIt->first] = 1;// set as an active nodes.

	 for (link_iterator linkIt = nod.inLinks.begin();
	 linkIt != nod.inLinks.end(); linkIt++)
	 isActives[linkIt->first] = 1;// set as an active nodes.
	 }
	 }
	 }
	 }
	 */

	vector<int>().swap(activeNodes);
	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	total_iterations_priorSPNodes++;

	delete[] randomGlobalArray;
	delete[] intSendPack;
	delete[] intAllRecvPack;
	delete[] doubleSendPack;
	delete[] doubleAllRecvPack;

//gettimeofday(&endPrio_moveSPnodes, NULL);

	/*total_time_prioritize_Spmove += elapsedTimeInSec(startPrio_moveSPnodes,
	 endPrio_moveSPnodes);*/

	endPrio_moveSPnodes = MPI_Wtime();
	total_time_prioritize_Spmove += endPrio_moveSPnodes - startPrio_moveSPnodes;

	return numMoved;
}

double Network::calculateModularityScore() {

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	map<int, int> module_indices;

	int index = 0;

	for (int i = 0; i < nNode; i++) {
		int modId = nodes[i].modIdx;
		if (module_indices.count(modId) > 0) {
			//do nothing
		} else {
			module_indices[modId] = index;
			index++;
		}
	}
	//printf("No. of module found in rank:%d is %d, total modules:%d\n", rank,
			//module_indices.size(), nModule);

	float* ai = new float[module_indices.size()]();
	float* eii = new float[module_indices.size()]();

	float Q = 0.0;

	for (int i = 0; i < nNode; i++) {

		Node& nd = nodes[i];
		int MdIndex = module_indices.at(nd.modIdx);

		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			int otherMdIndex = module_indices.at(otherNd.modIdx);

			if (nd.modIdx == otherNd.modIdx) {
				eii[otherMdIndex]++;
				ai[otherMdIndex]++;
			} else {
				ai[otherMdIndex]++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			if (nd.modIdx == otherNd.modIdx) {
				eii[MdIndex]++;
				ai[MdIndex]++;
			} else {
				ai[MdIndex]++;
			}
		}
	}

	for (int i = 0; i < module_indices.size(); i++) {
		eii[i] /= 2 * nEdge;
		ai[i] /= 2 * nEdge;
	}

	for (int i = 0; i < module_indices.size(); i++) {
		Q += eii[i] - pow(ai[i], 2);
	}

	delete[] eii;
	delete[] ai;

	return Q;

}

double Network::calculateConductance() {

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	map<int, int> module_indices;

	int index = 0;

	for (int i = 0; i < nNode; i++) {
		int modId = nodes[i].modIdx;
		if (module_indices.count(modId) > 0) {
			//do nothing
		} else {
			module_indices[modId] = index;
			index++;
		}
	}
	//printf("No. of module found in rank:%d is %d, total modules:%d\n", rank,
			//module_indices.size(), nModule);

	double intraEdge = 0, interEdge = 0;

	/*	float* ai = new float[module_indices.size()]();
	 float* eii = new float[module_indices.size()]();*/

	float C = 0.0;

	for (int i = 0; i < nNode; i++) {

		Node& nd = nodes[i];

		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			if (nd.modIdx == otherNd.modIdx) {
				intraEdge++;
			} else {
				interEdge++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			if (nd.modIdx == otherNd.modIdx) {
				intraEdge++;
			} else {
				interEdge++;
			}
		}
	}

	C = interEdge / (intraEdge + interEdge);

	return C;

}

double Network::calculateConductancePerModule() {

	int size;
	int rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ostringstream str1;
	str1 << rank;
	string a = "conductance_";
	string b = str1.str() + ".txt";
	string fileName = a + b;

	ofstream outputFile(fileName.c_str());

	map<int, int> module_indices;
	map<int, double> moduleWiseConductance;

	int index = 0;

	for (int i = 0; i < nNode; i++) {
		int modId = nodes[i].modIdx;
		if (module_indices.count(modId) > 0) {
			//do nothing
		} else {
			module_indices[modId] = index;
			index++;
		}
	}

	vector<int> intra(module_indices.size(), 0);
	vector<int> inter(module_indices.size(), 0);

	double C = 0.0;

	for (int i = 0; i < nNode; i++) {

		Node& nd = nodes[i];

		int mdIndex = module_indices.at(nd.modIdx);

		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			if (nd.modIdx == otherNd.modIdx) {
				intra[mdIndex]++;

			} else {
				inter[mdIndex]++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {

			Node& otherNd = nodes[linkIt->first];

			if (nd.modIdx == otherNd.modIdx) {
				intra[mdIndex]++;
			} else {
				inter[mdIndex]++;
			}
		}
	}

	for (int i = 0; i < module_indices.size(); i++) {
		if ((inter[i] + intra[i]) > 0) {
			moduleWiseConductance[i] = (double) inter[i]
					/ (inter[i] + intra[i]);

			/*			printf(
			 "pur rank:%d, i:%d, inter:%d, intra:%d, moduleWiseConductance:%f\n",
			 rank, i, inter[i], intra[i], moduleWiseConductance[i]);*/
		}
	}

	for (int i = 0; i < module_indices.size(); i++) {
		C += moduleWiseConductance[i];
		outputFile << i << " " << moduleWiseConductance[i] << endl;
		outputFile.flush();
	}

	outputFile.close();

	return C;

}

/*
 *	+ This is a parallel implementation of Network::moveSuperNodes() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */

int Network::parallelMoveSuperNodes(int numTh, double& tSequential) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nSuperNodes = superNodes.size();

// Generate random sequential order of nodes.
	vector<int> randomOrder(nSuperNodes);
	for (int i = 0; i < nSuperNodes; i++)
		randomOrder[i] = i;

	for (int i = 0; i < nSuperNodes; i++) {
		int target = R->randInt(nSuperNodes - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1; // This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);

	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nSuperNodes; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].sumPr > ndSize && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		if (bestResult.diffCodeLen < 0.0) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

					modules[newMod].numMembers += spMembers;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (ndDanglingSize > 0.0) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers -= spMembers;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved += spMembers;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

		}

	}	// END parallel for

	return numMoved;
}

/*
 *	+ This is a parallel implementation of Network::prioritize_moveSPnodes() function via OpenMP.
 *		We will loose the criteria for simple implementation and it may help to avoid local optima problem.
 */

int Network::prioritize_parallelMoveSPnodes(int numTh, double& tSequential,
		double vThresh) {

	struct timeval tStart, tEnd;

	gettimeofday(&tStart, NULL);

	int nActive = activeNodes.size();
	int nNextActive = 0;

	cout << "inside movesuper node" << endl;
// Generate random sequential order of nodes.
	vector<int> randomOrder(nActive);
	for (int i = 0; i < nActive; i++)
		randomOrder[i] = activeNodes[i];

	for (int i = 0; i < nActive; i++) {
		int target = R->randInt(nActive - 1);

		// swap numbers between i and target.
		int tmp = randomOrder[i];
		randomOrder[i] = randomOrder[target];
		randomOrder[target] = tmp;
	}

// Now randomOrder vector already had all activeNodes info,
// so we can reset activeNodes vector for adding for new active nodes.
	vector<int>().swap(activeNodes);

	int numMoved = 0;

	omp_set_num_threads(numTh);

	const int emptyTarget = nNode + 1;// This will be an indicator of moving to emptyModule.

	gettimeofday(&tEnd, NULL);

	tSequential += elapsedTimeInSec(tStart, tEnd);

#pragma omp parallel for
// Move each node to one of its neighbor modules in random sequential order.
	for (int i = 0; i < nActive; i++) {
		SuperNode& nd = superNodes[randomOrder[i]];	// look at i_th Node of the random sequential order.
		int oldMod = nd.ModIdx();

		unsigned int nModLinks = 0;	// The number of links to/from between the current node and other modules.

		flowmap outFlowToMod;	// <modID, flow> for outFlow...
		flowmap inFlowFromMod;	// <modID, flow> for inFlow...

		// count other modules that are connected to the current node.
		for (link_iterator linkIt = nd.outLinks.begin();
				linkIt != nd.outLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (outFlowToMod.count(newMod) > 0) {
				outFlowToMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = beta * linkIt->second;// initialization of the outFlow of the current newMod.
				inFlowFromMod[newMod] = 0.0;
				nModLinks++;
			}
		}

		for (link_iterator linkIt = nd.inLinks.begin();
				linkIt != nd.inLinks.end(); linkIt++) {
			int newMod = superNodes[linkIt->first].ModIdx();

			if (inFlowFromMod.count(newMod) > 0) {
				inFlowFromMod[newMod] += beta * linkIt->second;
			} else {
				outFlowToMod[newMod] = 0.0;
				inFlowFromMod[newMod] = beta * linkIt->second;
				nModLinks++;
			}
		}

		if (nModLinks != outFlowToMod.size())
			cout << "ALERT: nModLinks != outFlowToMod.size()." << endl;

		// copy node specific values for easy use and efficiency.
		double ndSize = nd.Size();	// p_nd.
		double ndTPWeight = nd.TeleportWeight();	// tau_nd.
		double ndDanglingSize = nd.DanglingSize();

		double oldExitPr1 = modules[oldMod].exitPr;
		double oldSumPr1 = modules[oldMod].sumPr;
		double oldSumDangling1 = modules[oldMod].sumDangling;
		double oldModTPWeight = modules[oldMod].sumTPWeight;

		double additionalTeleportOutFlow = (alpha * ndSize
				+ beta * ndDanglingSize) * (oldModTPWeight - ndTPWeight);
		double additionalTeleportInFlow = (alpha * (oldSumPr1 - ndSize)
				+ beta * (oldSumDangling1 - ndDanglingSize)) * ndTPWeight;

		// For teleportation and danling nodes.
		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			if (newMod == oldMod) {
				outFlowToMod[newMod] += additionalTeleportOutFlow;
				inFlowFromMod[newMod] += additionalTeleportInFlow;
			} else {
				outFlowToMod[newMod] += (alpha * ndSize + beta * ndDanglingSize)
						* modules[newMod].sumTPWeight;
				inFlowFromMod[newMod] += (alpha * modules[newMod].sumPr
						+ beta * modules[newMod].sumDangling) * ndTPWeight;
			}
		}

// Calculate flow to/from own module (default value if no links to own module).
		double outFlowToOldMod = additionalTeleportOutFlow;
		double inFlowFromOldMod = additionalTeleportInFlow;
		if (outFlowToMod.count(oldMod) > 0) {
			outFlowToOldMod = outFlowToMod[oldMod];
			inFlowFromOldMod = inFlowFromMod[oldMod];
		}

//////////////////// TODO: NEED TO IMPLEMENT THE OPTION TO MOVE TO EMPTY MODULE ////////////////
		if (modules[oldMod].sumPr > ndSize && emptyModules.size() > 0) {
			outFlowToMod[emptyTarget] = 0.0;
			inFlowFromMod[emptyTarget] = 0.0;
			nModLinks++;
		}

		MoveSummary currentResult;
		MoveSummary bestResult;

		double newExitPr1 = oldExitPr1 - nd.ExitPr() + outFlowToOldMod
				+ inFlowFromOldMod;

		bestResult.diffCodeLen = 0.0;// This is the default value, if we can't find diffCodeLen < 0, then don't move the node.

		for (flowmap::iterator it = outFlowToMod.begin();
				it != outFlowToMod.end(); it++) {
			int newMod = it->first;
			double outFlowToNewMod = it->second;
			double inFlowFromNewMod = inFlowFromMod[newMod];

			if (newMod != oldMod) {

// copy module specific values...
				double oldExitPr2 = modules[newMod].exitPr;
				double oldSumPr2 = modules[newMod].sumPr;

// Calculate status of current investigated movement of the node nd.
				currentResult.newModule = newMod;
				currentResult.sumPr1 = oldSumPr1 - ndSize;// This should be 0.0, because oldModule will be empty module.
				currentResult.sumPr2 = oldSumPr2 + ndSize;
				currentResult.exitPr1 = newExitPr1;
				currentResult.exitPr2 = oldExitPr2 + nd.ExitPr()
						- outFlowToNewMod - inFlowFromNewMod;

				currentResult.newSumExitPr = sumAllExitPr + newExitPr1
						+ currentResult.exitPr2 - oldExitPr1 - oldExitPr2;

// Calculate delta_L(M) = L(M)_new - L(M)_old
				double delta_allExit_log_allExit = pLogP(
						currentResult.newSumExitPr) - pLogP(sumAllExitPr);
				double delta_exit_log_exit = pLogP(currentResult.exitPr1)
						+ pLogP(currentResult.exitPr2) - pLogP(oldExitPr1)
						- pLogP(oldExitPr2);
				double delta_stay_log_stay = pLogP(
						currentResult.exitPr1 + currentResult.sumPr1)
						+ pLogP(currentResult.exitPr2 + currentResult.sumPr2)
						- pLogP(oldExitPr1 + oldSumPr1)
						- pLogP(oldExitPr2 + oldSumPr2);

// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
				currentResult.diffCodeLen = delta_allExit_log_allExit
						- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

				if (currentResult.diffCodeLen < bestResult.diffCodeLen) {
					// we need to update bestResult with currentResult.
					bestResult.diffCodeLen = currentResult.diffCodeLen;
					bestResult.newModule = currentResult.newModule;
				}
			}

		}

		// store the best possilbe movement information if necessary.
		//if (bestResult.diffCodeLen < 0.0) {
		if (bestResult.diffCodeLen < vThresh) {

			bool isEmptyTarget = false;
			bool validMove = true;// This will indicate the validity of the decided movement.
			int newMod = bestResult.newModule;
			int spMembers = nd.members.size();

#pragma omp critical (moveUpdate)
			{
				gettimeofday(&tStart, NULL);

// if newMod == emptyTarget, it indicates moves to empty module.
				if ((nEmptyMod > 0) && (newMod == emptyTarget)
						&& (modules[oldMod].numMembers > 1)) {
					newMod = emptyModules.back();
					isEmptyTarget = true;
				} else if (newMod == emptyTarget) {
					validMove = false;
				} else if (modules[newMod].numMembers == 0) {
					// This is the case that the algorithm thought there are some nodes in the new module since newMod != emptyTarget,
					// but the nodes are all moved to other modules so there are no nodes in there.
					// Thus, we don't know whether generating a new module will be better or not.
					// Therefore, don't move this option.
					//continue;
					validMove = false;
				}

				MoveSummary moveResult;

				if (validMove) {

					////////////////////////////////////////////////////////////////////
					/// THIS IS THE PART FOR EXAMINING QUALITY IMPROVEMENT OR NOT... ///
					////////////////////////////////////////////////////////////////////

					outFlowToOldMod = 0.0;
					double outFlowToNewMod = 0.0;
					inFlowFromOldMod = 0.0;
					double inFlowFromNewMod = 0.0;

					if (!isEmptyTarget) {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();

							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							} else if (toMod == newMod) {
								outFlowToNewMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();

							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							} else if (fromMod == newMod) {
								inFlowFromNewMod += beta * linkIt->second;
							}
						}
					} else {
						for (link_iterator linkIt = nd.outLinks.begin();
								linkIt != nd.outLinks.end(); linkIt++) {
							int toMod = superNodes[linkIt->first].ModIdx();
							if (toMod == oldMod) {
								outFlowToOldMod += beta * linkIt->second;
							}
						}

						for (link_iterator linkIt = nd.inLinks.begin();
								linkIt != nd.inLinks.end(); linkIt++) {
							int fromMod = superNodes[linkIt->first].ModIdx();
							if (fromMod == oldMod) {
								inFlowFromOldMod += beta * linkIt->second;
							}
						}
					}

					oldExitPr1 = modules[oldMod].exitPr;
					oldSumPr1 = modules[oldMod].sumPr;
					oldSumDangling1 = modules[oldMod].sumDangling;
					oldModTPWeight = modules[oldMod].sumTPWeight;

					// For teleportation and danling nodes.
					outFlowToOldMod += (alpha * ndSize + beta * ndDanglingSize)
							* (oldModTPWeight - ndTPWeight);
					inFlowFromOldMod += (alpha * (oldSumPr1 - ndSize)
							+ beta * (oldSumDangling1 - ndDanglingSize))
							* ndTPWeight;
					outFlowToNewMod += (alpha * ndSize + beta * ndDanglingSize)
							* modules[newMod].sumTPWeight;
					inFlowFromNewMod += (alpha * modules[newMod].sumPr
							+ beta * modules[newMod].sumDangling) * ndTPWeight;

					if (isEmptyTarget) {
						outFlowToNewMod = 0.0;
						inFlowFromNewMod = 0.0;
					}

					moveResult.exitPr1 = oldExitPr1 - nd.ExitPr()
							+ outFlowToOldMod + inFlowFromOldMod;

					// copy module specific values...
					double oldExitPr2 = modules[newMod].exitPr;
					double oldSumPr2 = modules[newMod].sumPr;

					// Calculate status of current investigated movement of the node nd.
					moveResult.newModule = newMod;
					moveResult.sumPr1 = oldSumPr1 - ndSize;	// This should be 0.0, because oldModule will be empty module.
					moveResult.sumPr2 = oldSumPr2 + ndSize;
					moveResult.exitPr2 = oldExitPr2 + nd.ExitPr()
							- outFlowToNewMod - inFlowFromNewMod;

					moveResult.newSumExitPr = sumAllExitPr + moveResult.exitPr1
							+ moveResult.exitPr2 - oldExitPr1 - oldExitPr2;

					// Calculate delta_L(M) = L(M)_new - L(M)_old
					double delta_allExit_log_allExit = pLogP(
							moveResult.newSumExitPr) - pLogP(sumAllExitPr);
					double delta_exit_log_exit = pLogP(moveResult.exitPr1)
							+ pLogP(moveResult.exitPr2) - pLogP(oldExitPr1)
							- pLogP(oldExitPr2);
					double delta_stay_log_stay = pLogP(
							moveResult.exitPr1 + moveResult.sumPr1)
							+ pLogP(moveResult.exitPr2 + moveResult.sumPr2)
							- pLogP(oldExitPr1 + oldSumPr1)
							- pLogP(oldExitPr2 + oldSumPr2);

					// delta_L(M) = delta_allExit - 2.0 * delta_exit_log_exit + delta_stay_log_stay.
					moveResult.diffCodeLen = delta_allExit_log_allExit
							- 2.0 * delta_exit_log_exit + delta_stay_log_stay;

					////////////////////////////////////
					////// THE END OF EXAMINATION //////
					////////////////////////////////////

					if (isEmptyTarget) {
						emptyModules.pop_back();

						nEmptyMod--;
						nModule++;
					}

					nd.setModIdx(newMod);
					for (int k = 0; k < spMembers; k++) {
						nd.members[k]->setModIdx(newMod);
					}

					modules[newMod].numMembers += spMembers;
					modules[newMod].exitPr = moveResult.exitPr2;
					modules[newMod].sumPr = moveResult.sumPr2;
					modules[newMod].stayPr = moveResult.exitPr2
							+ moveResult.sumPr2;
					modules[newMod].sumTPWeight += ndTPWeight;

					if (ndDanglingSize > 0.0) {
						modules[newMod].sumDangling += ndDanglingSize;
						modules[oldMod].sumDangling -= ndDanglingSize;
					}

					// update related to the oldMod...
					modules[oldMod].numMembers -= spMembers;
					modules[oldMod].exitPr = moveResult.exitPr1;
					modules[oldMod].sumPr = moveResult.sumPr1;
					modules[oldMod].stayPr = moveResult.exitPr1
							+ moveResult.sumPr1;
					modules[oldMod].sumTPWeight -= ndTPWeight;

					if (modules[oldMod].numMembers == 0) {
						nEmptyMod++;
						nModule--;
						emptyModules.push_back(oldMod);
					}

					sumAllExitPr = moveResult.newSumExitPr;
					codeLength += moveResult.diffCodeLen;

					numMoved += spMembers;
				}

				gettimeofday(&tEnd, NULL);

				tSequential += elapsedTimeInSec(tStart, tEnd);
			}	// END critical

// update activeNodes and isActives vectors.
// We have to add the following nodes in activeNodes: neighbors, members in oldMod & newMod.
// This can be done in parallel without any locking, since the written value is always same.
			for (link_iterator linkIt = nd.outLinks.begin();
					linkIt != nd.outLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;	// set as an active nodes.

			for (link_iterator linkIt = nd.inLinks.begin();
					linkIt != nd.inLinks.end(); linkIt++)
				isActives[linkIt->first] = 1;	// set as an active nodes.
		}

	}	// END parallel for

	for (int i = 0; i < isActives.size(); i++) {
		if (isActives[i] == 1) {
			activeNodes.push_back(i);
			isActives[i] = 0;	// reset the flag of isActives[i].
		}
	}

	return numMoved;
}

/**
 *	This function will update members vector in modules correspondingly.
 */
void Network::updateMembersInModule() 
{

	int size, rank;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto tik = std::chrono::high_resolution_clock::now();

	iterations_updateMembers++;

	int numModules = modules.size();

	vector<int>().swap(smActiveMods);
	smActiveMods.reserve(nModule);
	vector<int>().swap(lgActiveMods);
	lgActiveMods.reserve(nModule);

	for (int i = 0; i < numModules; i++) 
	{
		vector<Node*>().swap(modules[i].members);
		if (modules[i].numMembers > 10000)
			lgActiveMods.push_back(i);
		else if (modules[i].numMembers > 0)
			smActiveMods.push_back(i);
	}

	for (int i = 0; i < nNode; i++) 
	{
		modules[nodes[i].ModIdx()].members.push_back(&nodes[i]);
	}
	
	auto tok = std::chrono::high_resolution_clock::now();

	updateMembersTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();
}

void Network::updateSPMembersInModule() 
{
	int numModules = modules.size();

	for (int i = 0; i < numModules; i++)
		vector<Node*>().swap(modules[i].members);

	int numSPNodes = superNodes.size();
	for (int i = 0; i < numSPNodes; i++)
		modules[superNodes[i].ModIdx()].members.push_back(&superNodes[i]);
}

// calculate exit-probability based on node information and current module assignment.
void Network::updateCodeLength(int numTh, bool isSPNode) {

// calculate exit-probability for each module.
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

	omp_set_num_threads(numTh);

#pragma omp parallel for reduction (+:tempSumAllExit, exit_log_exit, stay_log_stay)
	for (unsigned int i = 0; i < nNode; i++) {
		if (modules[i].numMembers > 0) {
// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].sumTPWeight)
					* modules[i].sumPr;
			int curMod = modules[i].index;

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				Node * nd = modules[i].members[j];
				int nOutLinks = nd->outLinks.size();

				for (int k = 0; k < nOutLinks; k++) {
					if (!isSPNode
							&& nodes[nd->outLinks[k].first].ModIdx() != curMod)
						sumOutFlow += nd->outLinks[k].second;
					else if (isSPNode
							&& superNodes[nd->outLinks[k].first].ModIdx()
									!= curMod)
						sumOutFlow += nd->outLinks[k].second;
				}
			}

			exitPr += Network::beta
					* (sumOutFlow
							+ (1.0 - modules[i].sumTPWeight)
									* modules[i].sumDangling);

			modules[i].exitPr = exitPr;
			modules[i].stayPr = exitPr + modules[i].sumPr;

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(modules[i].stayPr);
		}
	}

	sumAllExitPr = tempSumAllExit;

	codeLength = pLogP(sumAllExitPr) - 2.0 * exit_log_exit + stay_log_stay
			- allNodes_log_allNodes;
}

// calculate exit-probability based on node information and current module assignment.
double Network::calculateCodeLength() {

	int size, rank;
	int startIndex, endIndex;
	const int arraySize = 3;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double* send_sum_exit_stay = new (nothrow) double[arraySize]();
//double receive_sum_exit_stay[arraySize] = { 0.0, 0.0, 0.0 };
//double* s_sum_exit_stay = new (nothrow) double[arraySize]();
	double* recv_sum_all_exit_stay = new (nothrow) double[arraySize * size]();

	if (!recv_sum_all_exit_stay || !send_sum_exit_stay) {
		printf(
				"ERROR in calculateCodeLength: process:%d does not have sufficient memory, process exiting...\n",
				rank);
		exit(1);
	}

	findAssignedPart(&startIndex, &endIndex, nNode, size, rank);

// calculate exit-probability for each module.
	double tempSumAllExit = 0.0;
	double exit_log_exit = 0.0;
	double stay_log_stay = 0.0;

	for (unsigned int i = startIndex; i < endIndex; i++) {
		if (modules[i].numMembers > 0) {
// exitPr w.r.t. teleportation.
			double exitPr = Network::alpha * (1.0 - modules[i].sumTPWeight)
					* modules[i].sumPr;
			int curMod = modules[i].index;

			double sumOutFlow = 0.0;
			int nMembers = modules[i].members.size();
			for (int j = 0; j < nMembers; j++) {
				int nOutLinks = modules[i].members[j]->outLinks.size();
				for (int k = 0; k < nOutLinks; k++) {
					if (nodes[modules[i].members[j]->outLinks[k].first].ModIdx()
							!= curMod)
						sumOutFlow += modules[i].members[j]->outLinks[k].second;
				}
			}

			exitPr += Network::beta
					* (sumOutFlow
							+ (1.0 - modules[i].sumTPWeight)
									* modules[i].sumDangling);

			tempSumAllExit += exitPr;
			exit_log_exit += pLogP(exitPr);
			stay_log_stay += pLogP(exitPr + modules[i].sumPr);
		}
	}

	send_sum_exit_stay[0] = tempSumAllExit;
	send_sum_exit_stay[1] = exit_log_exit;
	send_sum_exit_stay[2] = stay_log_stay;

	MPI_Allgather(send_sum_exit_stay, arraySize, MPI_DOUBLE,
			recv_sum_all_exit_stay, arraySize, MPI_DOUBLE,
			MPI_COMM_WORLD);

	for (int p = 0; p < size; p++) {
		if (rank != p) {
			int startIndex = p * arraySize;
			tempSumAllExit += recv_sum_all_exit_stay[startIndex];
			exit_log_exit += recv_sum_all_exit_stay[startIndex + 1];
			stay_log_stay += recv_sum_all_exit_stay[startIndex + 2];
		}
	}

	/*	for (int processId = 0; processId < size; processId++) {
	 if (processId != rank) {
	 MPI_Sendrecv(send_sum_exit_stay, 3, MPI_DOUBLE, processId, 0,
	 receive_sum_exit_stay, 3, MPI_DOUBLE, processId, 0,
	 MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

	 tempSumAllExit += receive_sum_exit_stay[0];
	 exit_log_exit += receive_sum_exit_stay[1];
	 stay_log_stay += receive_sum_exit_stay[2];
	 }
	 }*/

	double computedCodeLength = pLogP(tempSumAllExit) - 2.0 * exit_log_exit
			+ stay_log_stay - allNodes_log_allNodes;

	delete[] send_sum_exit_stay;
	delete[] recv_sum_all_exit_stay;
	return computedCodeLength;
}

// similar function of Greedy::level() in the original infomap_dir implementation.
// make a group of nodes as a Node (here, SuperNode)
// Granularity of the network is same but a group of nodes move together instead of moving individually.
void Network::convertModulesToSuperNodes(int numTh) {
//initialize superNodes vector for updating w.r.t. the current module status.

	auto tik = std::chrono::high_resolution_clock::now();

	auto ding = std::chrono::high_resolution_clock::now();

	iteration_convertModules++;

	vector<SuperNode>().swap(superNodes);
	superNodes.reserve(nModule);

	vector<unsigned int> modToSPNode(nNode);// indicate corresponding superNode ID from each module.

	int idx = 0;
	for (unsigned int i = 0; i < nNode; i++) 
	{
		if (modules[i].numMembers > 0) 
		{
			superNodes.push_back(SuperNode(modules[i], idx));
			modToSPNode[i] = idx;
			idx++;
		}
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();
	
	auto dong = std::chrono::high_resolution_clock::now();

	convertModulePart1 += std::chrono::duration_cast<std::chrono::nanoseconds>(dong - ding).count();

	auto ting = std::chrono::high_resolution_clock::now();

	omp_set_num_threads(numTh);


		/*
		int myID = omp_get_thread_num();
		int nTh = omp_get_num_threads();

		int start, end;

		findAssignedPart(&start, &end, numSPNode, nTh, myID);
		for (int i = start; i < end; i++) 
		{
		*/

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < numSPNode; i++) 
		{
			int numNodesInSPNode = superNodes[i].members.size();

			typedef unordered_map<int, double> EdgeMap;
			EdgeMap newOutLinks;

			for (int j = 0; j < numNodesInSPNode; j++) 
			{
				// Calculate newOutLinks from a superNode to other superNodes.
				Node* nd = superNodes[i].members[j];
				int nOutEdge = nd->outLinks.size();

				for (int k = 0; k < nOutEdge; k++) 
				{
					int toSPNode = modToSPNode[nodes[nd->outLinks[k].first].ModIdx()];

					if (toSPNode != i) 
					{	// outLinks to another superNode...
						pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(make_pair(toSPNode, nd->outLinks[k].second));
						if (!ret.second)
						{
							ret.first->second += nd->outLinks[k].second;
						}
					}
				}
			}

			superNodes[i].outLinks.reserve(newOutLinks.size());

			for (EdgeMap::iterator it = newOutLinks.begin(); it != newOutLinks.end(); it++) 
			{
				superNodes[i].outLinks.push_back(make_pair(it->first, it->second));
			}
		}

	auto tong = std::chrono::high_resolution_clock::now();

	convertModulePart2 += std::chrono::duration_cast<std::chrono::nanoseconds>(tong - ting).count();

	auto ping = std::chrono::high_resolution_clock::now();

// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) 
	{
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
		{
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(make_pair(i, superNodes[i].outLinks[j].second));
		}
	}

	auto pong = std::chrono::high_resolution_clock::now();

	convertModulePart3 += std::chrono::duration_cast<std::chrono::nanoseconds>(pong - ping).count();

	auto tok = std::chrono::high_resolution_clock::now();
	convert2SpNodeTime += std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();
}

void Network::displayOutlinksforSuperNodes() {

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numSPNode = superNodes.size();
	for (int it = 0; it < numSPNode; it++) {
		int noOutLinks = superNodes[it].outLinks.size();
		printf(
				"filtering by the process:%d, spNode:%d, number of outlinks:%d\n",
				rank, it, noOutLinks);
		for (int ite = 0; ite < noOutLinks; ite++) {
			printf(
					"filter by the process:%d, superNodes [%d].outLinks[%d]->otherEnd:%d, weight:%f\n",
					rank, it, ite, superNodes[it].outLinks[ite].first,
					superNodes[it].outLinks[ite].second);
		}
	}
}

void Network::displayInlinksforSuperNodes() {

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numSPNode = superNodes.size();
	for (int it = 0; it < numSPNode; it++) {
		int noinLinks = superNodes[it].inLinks.size();
		printf("refining by the process:%d, spNode:%d, number of inlinks:%d\n",
				rank, it, noinLinks);
		for (int ite = 0; ite < noinLinks; ite++) {
			printf(
					"refine by the process:%d, superNodes [%d].inLinks[%d]->otherEnd:%d, weight:%f\n",
					rank, it, ite, superNodes[it].inLinks[ite].first,
					superNodes[it].inLinks[ite].second);
		}
	}
}

//void Network::generateSuperNodesFromSubModules() {
void Network::generateSuperNodesFromSubModules(int numTh) {
//initialize superNodes vector for updating w.r.t. the current module status.
	vector<SuperNode>().swap(superNodes);

	int nSubMods = subModules.size();
	superNodes.reserve(nSubMods);

	for (int i = 0; i < nSubMods; i++) {
		superNodes.push_back(SuperNode(&subModules[i], i, *this));
	}

	/*
	 *	Calculate outLinks and inLinks between superNodes...
	 */
	int numSPNode = superNodes.size();

	typedef map<int, double> EdgeMap;

	omp_set_num_threads(numTh);

	for (int i = 0; i < numSPNode; i++) {
		SuperNode* spNode = &(superNodes[i]);

		int numNodesInSPNode = spNode->members.size();

		EdgeMap newOutLinks;

		for (int j = 0; j < numNodesInSPNode; j++) {
// Calculate newOutLinks from a superNode to other superNodes.
			Node* nd = spNode->members[j];
			int nOutEdge = nd->outLinks.size();

			for (int k = 0; k < nOutEdge; k++) {
				int toSPNode = ndToSubMod[nd->outLinks[k].first];

				if (toSPNode != i) {	// outLinks to another superNode...
					pair<EdgeMap::iterator, bool> ret = newOutLinks.insert(
							make_pair(toSPNode, nd->outLinks[k].second));
					if (!ret.second) {
						ret.first->second += nd->outLinks[k].second;
					}
				}
			}
		}

		spNode->outLinks.reserve(newOutLinks.size());

		double sumExitFlow = 0.0;

		for (EdgeMap::iterator it = newOutLinks.begin();
				it != newOutLinks.end(); it++) {
			spNode->outLinks.push_back(make_pair(it->first, it->second));
			sumExitFlow += it->second;
		}

		double teleportExitFlow = 0.0;
		for (int j = 0; j < numNodesInSPNode; j++) {
			Node* nd = spNode->members[j];
			teleportExitFlow += (alpha * nd->Size() + beta * nd->DanglingSize())
					* (1.0 - nd->TeleportWeight());
		}

		if (sumExitFlow == 0.0) {
			spNode->setExitPr(teleportExitFlow);
		} else {
			double exitProb = teleportExitFlow + beta * sumExitFlow;
			spNode->setExitPr(exitProb);
		}

	}

// update inLinks in SEQUENTIAL..
	for (int i = 0; i < numSPNode; i++) {
		int nOutLinks = superNodes[i].outLinks.size();
		for (int j = 0; j < nOutLinks; j++)
			superNodes[superNodes[i].outLinks[j].first].inLinks.push_back(
					make_pair(i, superNodes[i].outLinks[j].second));
	}
}

void Network::copyModule(Module * newM, Module * oldM) {
	newM->index = oldM->index;
	newM->exitPr = oldM->exitPr;
	newM->stayPr = oldM->stayPr;
	newM->sumPr = oldM->sumPr;
	newM->sumTPWeight = oldM->sumTPWeight;
	newM->sumDangling = oldM->sumDangling;
	newM->numMembers = oldM->numMembers;

	newM->members.clear();
	int nMembers = oldM->members.size();
	for (int i = 0; i < nMembers; i++)
		newM->members.push_back(oldM->members[i]);
}

void Network::copyModule(int newM, int oldM) {
	modules[newM].index = modules[oldM].index;
	modules[newM].exitPr = modules[oldM].exitPr;
	modules[newM].stayPr = modules[oldM].stayPr;
	modules[newM].sumPr = modules[oldM].sumPr;
	modules[newM].sumTPWeight = modules[oldM].sumTPWeight;
	modules[newM].sumDangling = modules[oldM].sumDangling;
	modules[newM].numMembers = modules[oldM].numMembers;

	modules[newM].members.clear();
	int nMembers = modules[oldM].members.size();
	for (int i = 0; i < nMembers; i++)
		modules[newM].members.push_back(modules[oldM].members[i]);
}

// Miscellaneous function..
void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID) {
	if (numNodes % numTh == 0) {
		*start = (numNodes / numTh) * myID;
		*end = (numNodes / numTh) * (myID + 1);
	} else {
		int block = numNodes / numTh;
		int modular = numNodes % numTh;

		if (myID < modular) {
			*start = (block + 1) * myID;
			*end = (block + 1) * (myID + 1);
		} else {
			*start = block * myID + modular;
			*end = block * (myID + 1) + modular;
		}
	}
}

