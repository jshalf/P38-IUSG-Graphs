#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <ctime>
#include <sys/time.h>
#include <omp.h>
#include "MersenneTwister.h"
#include "Node.h"
#include "Module.h"
#include "FileIO.h"
#include "timing.h"
#include <mpi.h>
#include </global/common/sw/cray/cnl7/haswell/metis/5.1.0/intel/19.0.3.199/nbtsmmb/include/metis.h>


using namespace std;

unsigned stou(char *s) {
	return strtoul(s, (char **) NULL, 10);
}

void stochastic_greedy_partition(Network &network, int numTh, double threshold,
		double vThresh, int maxIter, bool prior, bool fineTune, bool fast,
		bool inLoop, double& total_time_move, int& total_iterations_move,
		double& total_time_updateMembers, double& total_time_convertModules,
		double& total_time_prioritize_move,
		double& total_time_prioritize_Spmove, int& total_iterations_priorMove,
		int& total_iterations_priorMoveSP, int& total_iterations_updateMembers,
		int& total_iterations_convertModules, double& total_time_MPISendRecv,
		double& total_time_MPISendRecvSP,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules);
void partition_module_network(Network &network, int numTh, double threshold,
		int maxIter, bool fast, double& total_time_move,
		int& total_iterations_move, double& total_time_updateMembers,
		double& total_time_convertModules, int& total_iterations_updateMembers,
		int& total_iterations_convertModules,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules);
void generate_sub_modules(Network &network, int numTh, double threshold,
		int maxIter, double& total_time_move, int& total_iterations_move,
		double& total_time_updateMembers, double& total_time_convertModules,
		double& total_time_calibrate, int& total_iterations_updateMembers,
		int& total_iterations_convertModules,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules);
void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int>& origNodeID, int numTh, int iteration,
		double& total_time_calibrate);
void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int> &origNodeID, int iteration, double& total_time_calibrate);

void print_twoLevel_Cluster(Network network, string networkName, string outDir);

void findAssignedPart(int* start, int* end, int numNodes, int numTh, int myID);

int main(int argc, char *argv[]) {

	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 10) {
		cout
				<< "Call: ./ompRelaxmap <seed> <network.net> <# threads> <# attempts> <threshold> <vThresh> <maxIter> <outDir> <prior/normal>  [selflinks]"
				<< endl;
		exit(-1);
	}

	string outDir = string(argv[8]);
	int maxIter = atoi(argv[7]);	// Set the maximum number of iteration..
	int Ntrials = atoi(argv[4]);  // Set number of partition attempts
	int numThreads = atoi(argv[3]);	// Set number of threads...
	string line;
	string buf;

	MTRand *R = new MTRand(stou(argv[1]));

	string infile = string(argv[2]);
	string networkFile = string(argv[2]);
	string networkName(networkFile.begin() + networkFile.find_last_of("/"),
			networkFile.begin() + networkFile.find_last_of("."));
	string networkType(infile.begin() + infile.find_last_of("."), infile.end());

	double threshold = atof(argv[5]);
	double vThresh = atof(argv[6]);	// vertex-threshold: threshold for each vertex-movement.

	cout << "Threshold = " << threshold << ", Vertex-Threshold = " << vThresh
			<< endl;
	vThresh *= -1;	// change the threshold value for negative.

	string priorFlag = string(argv[9]);

	bool prior = false;
	if (priorFlag == "prior")
		prior = true;

	string metisFile = string(argv[10]);

	bool includeSelfLinks = false;
	if (argc == 12) {
		string selfLinks(argv[11]);
		if (selfLinks == "selflinks")
			includeSelfLinks = true;
	}

	Network origNetwork;	// default constructor is called.

	origNetwork.R = R;

	// time values for measuring elapsed times for each step...
	struct timeval allStart, allEnd;
	struct timeval noIOstart, noIOend;
	struct timeval start, end;

	gettimeofday(&allStart, NULL);
	gettimeofday(&start, NULL);

	if (networkType == ".net") {
		load_pajek_format_network(networkFile, origNetwork);
	} else if (networkType == ".csr") {
		load_csr_format_network(networkFile, origNetwork, metisFile);
	} else {
		load_linkList_format_network(networkFile, origNetwork);
	}
	//MPI_Send(&origNetwork,1,)
	gettimeofday(&end, NULL);

	cout << "Time for reading input data : " << elapsedTimeInSec(start, end)
			<< " (sec)" << endl;

	gettimeofday(&noIOstart, NULL);

	int nNode = origNetwork.NNode();
	double totNodeWeights = origNetwork.TotNodeWeights();

	cout << "total Node weights = " << totNodeWeights << endl;

	//here let's build the parmetis graph
	//origNetwork.buildParMetis();

	gettimeofday(&start, NULL);

	for (int i = 0; i < nNode; i++) {
		origNetwork.nodes[i].setTeleportWeight(
				origNetwork.nodes[i].NodeWeight() / totNodeWeights);
	}

	int NselfLinks = 0;
	for (map<pair<int, int>, double>::iterator it = origNetwork.Edges.begin();
			it != origNetwork.Edges.end(); it++) {

		int from = it->first.first;
		int to = it->first.second;
		double weight = it->second;
		if (weight > 0.0) {
			if (from == to) {
				NselfLinks++;
				//if(includeSelfLinks)
				//	origNetwork.nodes[from]->selfLink += weight;
			} else {
				origNetwork.nodes[from].outLinks.push_back(
						make_pair(to, weight));
				// we will going to update inLinks, after we got final flow of the network.
				//origNetwork.nodes[to].inLinks.push_back(make_pair(from,weight));
			}
		}
	}

	if (includeSelfLinks)
	//cout << "including " <<  NselfLinks << " self link(s)." << endl;
	{
		cout << "current version always excludes self links.\nignoring "
				<< NselfLinks << " self link(s)." << endl;
	} else {
		cout << "ignoring " << NselfLinks << " self link(s)." << endl;
	}

	//Swap vector to free memory
	map<pair<int, int>, double>().swap(origNetwork.Edges);

	cout << "DONE: Parsing the given network  ..." << endl;

	gettimeofday(&end, NULL);
	cout << "Time for parsing the given network : "
			<< elapsedTimeInSec(start, end) << " (sec)" << endl;

	gettimeofday(&start, NULL);
	// Initialization.

	double total_time_initiate = 0.0;
	double total_time_calibrate = 0.0;
	double total_time_updateMembers = 0.0;
	double total_time_convertModules = 0.0;
	double total_time_move = 0.0;
	double total_time_calcCodelen = 0.0;
	double total_time_prioritize_move = 0.0;
	double total_time_prioritize_Spmove = 0.0;
	double total_time_MPISendRecv = 0.0;
	double total_time_MPISendRecvSP = 0.0;
	double total_time_MPISendRecvUpdateMem = 0.0;
	double total_time_MPISendRecvConvertModules = 0.0;
	int total_iterations_move = 0;
	int total_iterations_priorMove = 0;
	int total_iterations_priorMoveSP = 0;
	int total_iterations_updateMembers = 0;
	int total_iterations_convertModules = 0;

	origNetwork.initiate(numThreads, total_time_initiate, total_time_calibrate);

	// Now update inLinks..
	for (int i = 0; i < nNode; i++) {
		int nOutLinks = origNetwork.nodes[i].outLinks.size();

		for (int j = 0; j < nOutLinks; j++)
			origNetwork.nodes[origNetwork.nodes[i].outLinks[j].first].inLinks.push_back(
					make_pair(i, origNetwork.nodes[i].outLinks[j].second));

	}

	gettimeofday(&end, NULL);

	printf("DONE: Initiate() ... for rank:%d in %f sec\n", rank,
			elapsedTimeInSec(start, end));
	cout << "Initial Code Length: " << origNetwork.CodeLength() / log(2.0)
			<< " in " << origNetwork.NModule() << " modules.\n";

	// copy size of each node for print in order.
	vector<double> nodeSize(nNode);
	for (int i = 0; i < nNode; i++) {
		nodeSize[i] = origNetwork.nodes[i].Size();
	}
	cout << "Now partition the network starts...\n";
	gettimeofday(&start, NULL);

	bool fineTune = true;
	bool fast = false;	// This will be true only for sub-module partitioning...

	int step = 1;

	// Initial SuperStep running ...
	double oldCodeLength = origNetwork.CodeLength();

	stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh,
			maxIter, prior, fineTune, fast, false, total_time_move,
			total_iterations_move, total_time_updateMembers,
			total_time_convertModules, total_time_prioritize_move,
			total_time_prioritize_Spmove, total_iterations_priorMove,
			total_iterations_priorMoveSP, total_iterations_updateMembers,
			total_iterations_convertModules, total_time_MPISendRecv,
			total_time_MPISendRecvSP, total_time_MPISendRecvUpdateMem,
			total_time_MPISendRecvConvertModules);

	/*	cout << "SuperStep [" << step << "] - codeLength = "
	 << origNetwork.CodeLength() / log(2.0) << " in "
	 << origNetwork.NModule() << " modules." << endl;*/

	printf("SuperStep[%d] - codeLength = %f in %d modules.\n", step,
			origNetwork.CodeLength() / log(2.0), origNetwork.NModule());

	bool nextIter = true;

	if ((oldCodeLength - origNetwork.CodeLength()) / log(2.0) < threshold
			|| origNetwork.CodeLength() <= 0.0) {
		nextIter = false;
		/*		printf("1new network for rank:%d with codelength:%f\n", rank,
		 origNetwork.CodeLength());*/
	}

	struct timeval subStart, subEnd;

	//nextIter = false; // I am forcefully doing it to stop executing the while loop below for now

	while (nextIter) {

		oldCodeLength = origNetwork.CodeLength();

		stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh,
				maxIter, prior, fineTune, fast, true, total_time_move,
				total_iterations_move, total_time_updateMembers,
				total_time_convertModules, total_time_prioritize_move,
				total_time_prioritize_Spmove, total_iterations_priorMove,
				total_iterations_priorMoveSP, total_iterations_updateMembers,
				total_iterations_convertModules, total_time_MPISendRecv,
				total_time_MPISendRecvSP, total_time_MPISendRecvUpdateMem,
				total_time_MPISendRecvConvertModules);

		step++;

		printf("SuperStep[%d] - codeLength = %f in %d modules.\n", step,
				origNetwork.CodeLength() / log(2.0), origNetwork.NModule());

		/*		cout << "SuperStep [" << step << "] - codeLength = "
		 << origNetwork.CodeLength() / log(2.0) << " in "
		 << origNetwork.NModule() << " modules." << endl;*/

		if ((oldCodeLength - origNetwork.CodeLength()) / log(2.0) < threshold
				|| origNetwork.CodeLength() < 0.0) {
			nextIter = false;
			/*			printf("2new network for rank:%d with codelength:%f\n", rank,
			 origNetwork.CodeLength());*/
		}

		//TODO: for now I am ignoring the operation of reconstructing network based on found modules

		/*fineTune = !fineTune; // fine-tune and coarse-tune will be done alternatively.

		 if (nextIter && !fineTune) {
		 // Next iteration will be Coarse Tune.
		 gettimeofday(&subStart, NULL);

		 generate_sub_modules(origNetwork, numThreads, threshold, maxIter,
		 total_time_move, total_iterations_move,
		 total_time_updateMembers, total_time_convertModules,
		 total_time_calibrate, total_iterations_updateMembers,
		 total_iterations_convertModules);

		 gettimeofday(&subEnd, NULL);
		 printf("Time for finding sub-modules for rank:%d is:%f sec\n", rank,
		 elapsedTimeInSec(subStart, subEnd));
		 }*/
	}

	gettimeofday(&end, NULL);
	cout << "Time for partitioning : " << elapsedTimeInSec(start, end)
			<< " (sec)" << endl;

	cout << "DONE: Code Length = " << origNetwork.CodeLength() / log(2.0)
			<< " in ";
	cout << origNetwork.NModule() << " modules, with " << nNode << " nodes.\n"
			<< endl;

	gettimeofday(&noIOend, NULL);
	gettimeofday(&allEnd, NULL);

	cout << "Overall Elapsed Time for Module Detection (w/o file IO): "
			<< elapsedTimeInSec(noIOstart, noIOend) << " (sec)" << endl;
	cout << "Overall Elapsed Time for Module Detection (w/ file Reading): "
			<< elapsedTimeInSec(allStart, allEnd) << " (sec)" << endl;

	cout << "\nComputed Code Length = "
			<< origNetwork.calculateCodeLength(total_time_calcCodelen)
					/ log(2.0) << endl;

/*	double modularity = origNetwork.calculateModularityScore();

	printf("\nmodularity score for the network:%f\n", modularity);

	double conductance = origNetwork.calculateConductance();

	printf("\nconductance of the network:%f\n\n", conductance);*/

	/*	double cumulativeConductance = origNetwork.calculateConductancePerModule();

	 printf("\ncumulativeConductance of the network:%f\n\n",
	 cumulativeConductance);*/

	double communication_time = total_time_MPISendRecv
			+ total_time_MPISendRecvSP;

	double total_time = total_time_prioritize_move
			+ total_time_prioritize_Spmove;

	double computation_time = total_time - communication_time;

	printf("time for initiate function in rank:%d is %f\n", rank,
			total_time_initiate);

	printf("time for calibrate in rank:%d is %f\n", rank, total_time_calibrate);

	printf(
			"time for updateMembersInModule in rank:%d is %f for total:%d iterations\n",
			rank, total_time_updateMembers, total_iterations_updateMembers);

	printf(
			"time for convertModulesToSuperNodes in rank:%d is %f for total:%d iterations\n",
			rank, total_time_convertModules, total_iterations_convertModules);

	printf("time for calculateCodeLength in rank:%d is %f\n", rank,
			total_time_calcCodelen);

	printf("time for move in rank:%d is %f for total:%d iterations\n", rank,
			total_time_move, total_iterations_move);

	printf(
			"time for prioritize_move in rank:%d is %f for total:%d iterations\n",
			rank, total_time_prioritize_move, total_iterations_priorMove);

	printf(
			"time for MPI_SendRecv in prioritize_move in rank:%d is %f for total:%d iterations\n",
			rank, total_time_MPISendRecv, total_iterations_priorMove);

	printf(
			"time for prioritize_Spmove in rank:%d is %f for total:%d iterations\n",
			rank, total_time_prioritize_Spmove, total_iterations_priorMoveSP);

	printf(
			"time for MPI_SendRecv in prioritize_Spmove in rank:%d is %f for total:%d iterations\n",
			rank, total_time_MPISendRecvSP, total_iterations_priorMoveSP);

	printf("time for computataion in rank:%d is %f\n", rank, computation_time);
	printf("time for MPI_communication in rank:%d is %f\n", rank,
			communication_time);

//Print two-level clustering result in .tree file
	print_twoLevel_Cluster(origNetwork, networkName, outDir);

// Print partition in Pajek's .clu format
	ofstream outFile;
	ostringstream oss;

	oss.str("");
	oss << outDir << "/" << networkName << ".clu";
	outFile.open(oss.str().c_str());
	outFile << "*Vertices " << nNode << "\x0D\x0A";
	for (int i = 0; i < nNode; i++) {
		outFile << origNetwork.nodes[i].ModIdx() + 1 << "\x0D\x0A";
	}
	outFile.close();

	MPI_Finalize();
}

/*
 * Procedure will be following:
 *	1) in random sequential order, each node is moved to its neighbor module that results in the largest gain of the map eq.
 *	   If no move results in a gain of the map equation, the node stays in its original module.
 *	2) repeated 1) procedure, each time in a new random sequential order, until no move generates a gain of the map EQ.
 *
 *	The 1) procedure is implemented in Network::move() function.
 */
void stochastic_greedy_partition(Network &network, int numTh, double threshold,
		double vThresh, int maxIter, bool prior, bool fineTune, bool fast,
		bool inLoop, double& total_time_move, int& total_iterations_move,
		double& total_time_updateMembers, double& total_time_convertModules,
		double& total_time_prioritize_move,
		double& total_time_prioritize_Spmove, int& total_iterations_priorMove,
		int& total_iterations_priorMoveSP, int& total_iterations_updateMembers,
		int& total_iterations_convertModules, double& total_time_MPISendRecv,
		double& total_time_MPISendRecvSP,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules) {

	double oldCodeLength = network.CodeLength();
	int iter = 0;
	bool stop = false;
	int rank, size;
	int tag = 0; //this tag is to track the call of prioritize_moveSPNode from while loop or do-while loop

	int first = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	struct timeval outer_T1, outer_T2;
	struct timeval inner_T1, inner_T2;
	struct timeval seq_T1, seq_T2;
	struct timeval convert_T1, convert_T2;

	double tSequential = 0.0;

	gettimeofday(&outer_T1, NULL);

	int nActiveUnits = (fineTune) ? network.NNode() : network.superNodes.size();

// set initial active nodes list ...
	vector<char>(nActiveUnits).swap(network.isActives);
	vector<int>(nActiveUnits).swap(network.activeNodes);
	map<int, int>().swap(network.activeVertices);

	for (int i = 0; i < nActiveUnits; i++) {
		network.activeNodes[i] = i;
		network.isActives[i] = 0;	// initially set inactive nodes.
		network.activeVertices.insert(make_pair(i, 1));
	}

	int numMoved = 0;

	while (!stop && iter < maxIter) {
		gettimeofday(&inner_T1, NULL);

		oldCodeLength = network.CodeLength();

		if (fineTune) {
			//if (numTh == 1) {
			if (prior) {
				numMoved = network.prioritize_move(vThresh, iter, inLoop,
						total_time_prioritize_move, total_iterations_priorMove,
						total_time_MPISendRecv);

			} else {
				numMoved = network.move(iter, total_time_move,
						total_iterations_move);
			}
			/*} else {
			 if (prior) {
			 numMoved = network.prioritize_parallelMove(numTh,
			 tSequential, vThresh);
			 } else {
			 numMoved = network.parallelMove(numTh, tSequential);
			 }
			 }*/
		} else {
			//if (numTh == 1) {
			if (prior) {
				numMoved = network.prioritize_moveSPnodes(vThresh, numTh, iter,
						inLoop, total_time_prioritize_Spmove,
						total_iterations_priorMoveSP, total_time_MPISendRecvSP);
			} else {
				numMoved = network.moveSuperNodes(iter);// If at least one node is moved, return true. Otherwise, return false.
			}
			/*} else {
			 if (prior) {
			 numMoved = network.prioritize_parallelMoveSPnodes(numTh,
			 tSequential, vThresh);
			 } else {
			 numMoved = network.parallelMoveSuperNodes(numTh,
			 tSequential);
			 }
			 }*/
		}
		iter++;
		/*
		 if (inLoop == false) {
		 printf(
		 "the old codelength for rank:%d, is:%f, the new codeLength is:%f, the iteration:%d blah blah blah\n",
		 rank, oldCodeLength, network.CodeLength() / log(2.0), iter);
		 }*/

		if (network.CodeLength() <= 0
				|| ((oldCodeLength - network.CodeLength()) / log(2.0)
						< threshold)) {
			stop = true;	//moved = false;
		}
		gettimeofday(&inner_T2, NULL);

		/*		// Print code length per iteration for DEBUG purpose.
		 cout << "Phor rank:" << rank << " Iteration " << iter
		 << ": code length =\t" << network.CodeLength() / log(2.0)
		 << "\t, ";
		 cout << "elapsed time:\t" << elapsedTimeInSec(inner_T1, inner_T2)
		 << "\t(sec), ";
		 cout << "accumulated time:\t" << elapsedTimeInSec(outer_T1, inner_T2)
		 << "\t(sec)\t";
		 cout << "sumExitPr = " << network.SumAllExitPr() << "\t";
		 cout << "numMoved:\t" << numMoved << endl;*/
		/*		cout << "Iteration " << iter << ": code length =\t"
		 << network.CodeLength() / log(2.0) << "\t, ";*/

		//double codel = network.CodeLength() / log(2.0);
		double unnecessary_time = 0.0;

		/*printf(
				"For rank:%d, Iteration:%d, code length:%f, numMoved:%d, number of active nodes:%ld\n",
				rank, iter,
				network.calculateCodeLength(unnecessary_time) / log(2.0),
				numMoved, network.activeNodes.size());*/
	}

	int outerLoop = 1;

	network.updateMembersInModule(total_time_updateMembers,
			total_iterations_updateMembers);

	tSequential += elapsedTimeInSec(seq_T1, seq_T2);

	if (fast)
		return;

	double tConvert = 0.0;

	do {

		oldCodeLength = network.CodeLength();
		stop = false;

		network.convertModulesToSuperNodes(numTh, total_time_convertModules,
				total_iterations_convertModules,
				total_time_MPISendRecvConvertModules);

		tConvert += elapsedTimeInSec(convert_T1, convert_T2);

		nActiveUnits = network.superNodes.size();

		// set initial active nodes list ...
		vector<char>(nActiveUnits).swap(network.isActives);
		vector<int>(nActiveUnits).swap(network.activeNodes);
		for (int i = 0; i < nActiveUnits; i++) {
			network.activeNodes[i] = i;	// initially all active units are active.
			network.isActives[i] = 0;	// initially set inactive nodes.
		}

		int spIter = 0;
		while (!stop && spIter < maxIter) {
			gettimeofday(&inner_T1, NULL);

			double innerOldCodeLength = network.CodeLength();

			//if (numTh == 1) {
			if (prior) {
				numMoved = network.prioritize_moveSPnodes(vThresh, numTh, spIter,
						inLoop, total_time_prioritize_Spmove,
						total_iterations_priorMoveSP, total_time_MPISendRecvSP);
			} else
				numMoved = network.moveSuperNodes(spIter);
			/*} else {
			 if (prior)
			 numMoved = network.prioritize_parallelMoveSPnodes(numTh,
			 tSequential, vThresh);
			 else
			 numMoved = network.parallelMoveSuperNodes(numTh,
			 tSequential);
			 }*/

			spIter++;

			if (inLoop == false) {
				/*				printf(
				 "the old codelength for rank:%d, in superNode is:%f, the new codeLength in superNode is:%f, the iteration:%d\n",
				 rank, innerOldCodeLength, network.CodeLength(), spIter);*/
			}
			if (network.CodeLength() <= 0.0
					|| (innerOldCodeLength - network.CodeLength()) / log(2.0)
							< threshold) {
				stop = true;	//moved = false;
			}
			gettimeofday(&inner_T2, NULL);

			/*			// Print code length per spIter for DEBUG purpose.
			 cout << "Phor rank:" << rank << " SuperIteration " << outerLoop
			 << "-" << spIter << ": code length =\t"
			 << network.CodeLength() / log(2.0) << "\t, ";
			 cout << "elapsed time:\t" << elapsedTimeInSec(inner_T1, inner_T2)
			 << "\t(sec), ";
			 cout << "accumulated time:\t"
			 << elapsedTimeInSec(outer_T1, inner_T2) << "\t(sec)\t";
			 cout << "sumExitPr = " << network.SumAllExitPr() << "\t";
			 cout << "numMoved:\t" << numMoved << endl;*/

			double unnecessary_time = 0.0;

/*			printf(
					"For rank:%d, SuperIteration:%d, code length:%f, numMoved:%d\n",
					rank, spIter,
					network.calculateCodeLength(unnecessary_time) / log(2.0),
					numMoved);*/
		}

		gettimeofday(&seq_T1, NULL);

		network.updateMembersInModule(total_time_updateMembers,
				total_iterations_updateMembers);

		outerLoop++;

		tag++;

		gettimeofday(&seq_T2, NULL);
		tSequential += elapsedTimeInSec(seq_T1, seq_T2);

	} while ((oldCodeLength - network.CodeLength()) / log(2.0) > threshold);

	gettimeofday(&outer_T2, NULL);

	/*	cout << "Sequential running time for partition: " << tSequential << " (sec)"
	 << endl;
	 cout << "Time for converting Module to SuperNode: " << tConvert << " (sec)"
	 << endl;
	 cout << "Overall time for partition: "
	 << elapsedTimeInSec(outer_T1, outer_T2) << "\t(sec)" << endl;*/
}

/**
 *	This function will be called for partitioning sub-Module of each module of the original graph.
 *	Thus, we would like to reduce printing from this function for providing high-level log.
 */
void partition_module_network(Network &network, int numTh, double threshold,
		int maxIter, bool fast, double& total_time_move,
		int& total_iterations_move, double& total_time_updateMembers,
		double& total_time_convertModules, int& total_iterations_updateMembers,
		int& total_iterations_convertModules,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules) {

	int rank, size;

	struct timeval startPartition, endPartition;

	gettimeofday(&startPartition, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double oldCodeLength = network.CodeLength();

	int iter = 0;
	bool stop = false;
	double tSequential;

	int numMoved = 0;

	while (!stop && iter < maxIter) {
		oldCodeLength = network.CodeLength();

		if (numTh == 1) {
			numMoved = network.move(iter, total_time_move,
					total_iterations_move);
		} else {
			numMoved = network.parallelMove(numTh, tSequential);
		}

		iter++;

		if ((oldCodeLength - network.CodeLength()) / log(2.0) < threshold)
			stop = true;
	}

	int outerLoop = 1;

	network.updateMembersInModule(total_time_updateMembers,
			total_iterations_updateMembers);

	if (fast)
		return;

	do {
		oldCodeLength = network.CodeLength();

		stop = false;
		network.convertModulesToSuperNodes(numTh, total_time_convertModules,
				total_iterations_convertModules,
				total_time_MPISendRecvConvertModules);

		int spIter = 0;
		while (!stop && spIter < maxIter) {
			double innerOldCodeLength = network.CodeLength();

			if (numTh == 1) {
				numMoved = network.moveSuperNodes(spIter);
			} else {
				numMoved = network.parallelMoveSuperNodes(numTh, tSequential);
			}

			spIter++;

			if ((innerOldCodeLength - network.CodeLength()) / log(2.0)
					< threshold)
				stop = true;
		}

		network.updateMembersInModule(total_time_updateMembers,
				total_iterations_updateMembers);

		outerLoop++;

	} while ((oldCodeLength - network.CodeLength()) / log(2.0) > threshold);

	gettimeofday(&endPartition, NULL);

	//printf("time for partition_module_network in rank:%d is %f\n", rank,

}

void generate_sub_modules(Network &network, int numTh, double threshold,
		int maxIter, double& total_time_move, int& total_iterations_move,
		double& total_time_updateMembers, double& total_time_convertModules,
		double& total_time_calibrate, int& total_iterations_updateMembers,
		int& total_iterations_convertModules,
		double& total_time_MPISendRecvUpdateMem,
		double& total_time_MPISendRecvConvertModules) {

	int rank, size;

	struct timeval startGenModules, endGenModules;

	gettimeofday(&startGenModules, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numNodes = network.NNode();

	struct timeval t1, t2;

	gettimeofday(&t1, NULL);

	vector<SubModule>().swap(network.subModules);// swap network.subModules vector with the empty vector.
	network.subModules.reserve(numNodes);

	vector<int>(numNodes).swap(network.ndToSubMod);

	gettimeofday(&t2, NULL);
	cout << "initialization time for generate_sub_modules(): "
			<< elapsedTimeInSec(t1, t2) << endl;

	struct timeval tPar1, tPar2;

	gettimeofday(&tPar1, NULL);

	MTRand *Rand = new MTRand();

	int numSmallMods = network.smActiveMods.size();
	cout << "number of small modules: " << numSmallMods << endl;

	for (int i = 0; i < numSmallMods; i++) {

		int iteration = 0;

		Module* mod = &(network.modules[network.smActiveMods[i]]);
		// check whether the current module has more than one node or not.
		if (mod->numMembers > 1) {
			int modIdx = mod->index;

			map<int, int> origNodeID;//map from newNodeID to original Node ID. a.k.a. <newNodeID, origNodeID>

			Network newNetwork;

			generate_network_from_module(newNetwork, mod, origNodeID, iteration,
					total_time_calibrate);
			newNetwork.R = Rand;

			partition_module_network(newNetwork, 1, threshold, maxIter, true,
					total_time_move, total_iterations_move,
					total_time_updateMembers, total_time_convertModules,
					total_iterations_updateMembers,
					total_iterations_convertModules,
					total_time_MPISendRecvUpdateMem,
					total_time_MPISendRecvConvertModules);	// fast = true..

			int nActiveMods = newNetwork.smActiveMods.size();
			// Adding sub-modules from a new network of the corresponding module to the list of subModules...
			for (int j = 0; j < nActiveMods; j++) {
				SubModule subMod(newNetwork.modules[newNetwork.smActiveMods[j]],
						origNodeID, modIdx);
				network.subModules.push_back(subMod);
			}
		} else {
			// This is the special case that the module has ONLY ONE member.
			SubModule subMod(*mod);
			network.subModules.push_back(subMod);
		}

		iteration++;

	}

	gettimeofday(&tPar2, NULL);

	cout << "Time for parallel for loop for SMALL-MODULES:\t"
			<< elapsedTimeInSec(tPar1, tPar2) << " (sec)" << endl;

///////////////////////////////
// Larger-Modules

	gettimeofday(&tPar1, NULL);

	int numLargeMods = network.lgActiveMods.size();
	cout << "number of large modules: " << numLargeMods << endl;

	for (int i = 0; i < numLargeMods; i++) {
		Module* mod = &(network.modules[network.lgActiveMods[i]]);

		int iteration = 0;

		// NO-NEED to check the size of the current module.
		int modIdx = mod->index;

		map<int, int> origNodeID; //map from newNodeID to original Node ID. a.k.a. <newNodeID, origNodeID>

		Network newNetwork;
		generate_network_from_module(newNetwork, mod, origNodeID, numTh,
				iteration, total_time_calibrate);
		newNetwork.R = Rand;

		partition_module_network(newNetwork, numTh, threshold, maxIter, true,
				total_time_move, total_iterations_move,
				total_time_updateMembers, total_time_convertModules,
				total_iterations_updateMembers, total_iterations_convertModules,
				total_time_MPISendRecvUpdateMem,
				total_time_MPISendRecvConvertModules); // fast = true..

		// Adding sub-modules from a new network of the corresponding module to the list of subModules...
		int nActiveMods = newNetwork.smActiveMods.size();

		for (int j = 0; j < nActiveMods; j++) {
			SubModule subMod(newNetwork.modules[newNetwork.smActiveMods[j]],
					origNodeID, modIdx);
			network.subModules.push_back(subMod);
		}

		nActiveMods = newNetwork.lgActiveMods.size();

		for (int j = 0; j < nActiveMods; j++) {
			SubModule subMod(newNetwork.modules[newNetwork.lgActiveMods[j]],
					origNodeID, modIdx);
			network.subModules.push_back(subMod);
		}

		iteration++;
	}

	gettimeofday(&tPar2, NULL);

	cout << "Time for parallel for loop for LARGE-MODULES:\t"
			<< elapsedTimeInSec(tPar1, tPar2) << " (sec)" << endl;

	gettimeofday(&t1, NULL);

	int numSubMods = 0;

	for (vector<SubModule>::iterator it = network.subModules.begin();
			it != network.subModules.end(); it++) {
		for (vector<int>::iterator ndIt = it->members.begin();
				ndIt != it->members.end(); ndIt++) {
			network.ndToSubMod[*ndIt] = numSubMods;
		}
		numSubMods++;
	}

	gettimeofday(&t2, NULL);
	cout << "sequential subModules push_back() time:\t"
			<< elapsedTimeInSec(t1, t2) << " (sec)" << endl;

	gettimeofday(&t1, NULL);
	network.generateSuperNodesFromSubModules(numTh);
	gettimeofday(&t2, NULL);

	gettimeofday(&endGenModules, NULL);

	printf("time for generate_sub_modules in rank:%d is %f\n", rank,
			elapsedTimeInSec(startGenModules, endGenModules));

	cout << "generateSuperNodesFromSubModules() time:\t"
			<< elapsedTimeInSec(t1, t2) << " (sec)" << endl;
}

void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int> &origNodeID, int iteration,
		double& total_time_calibrate) {

	int rank, size;

	struct timeval startGenNetworkmodule, endGenNetworkmodule;

	gettimeofday(&startGenNetworkmodule, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numMembers = mod->numMembers;

	newNetwork.modules = vector<Module>(numMembers);

	map<int, int> newNodeID;	// key = origNodeID --> value =  newNodeID.

	int newIdx = 0;
	for (vector<Node *>::iterator it = mod->members.begin();
			it != mod->members.end(); it++) {
		newNodeID[(*it)->ID()] = newIdx;
		origNodeID[newIdx] = (*it)->ID();
		Node nd(newIdx, (*it)->Size());
		nd.setNodeWeight((*it)->NodeWeight());
		nd.setTeleportWeight((*it)->TeleportWeight());
		nd.setDanglingSize((*it)->DanglingSize());
		nd.setIsDangling((*it)->IsDangling());
		newNetwork.nodes.push_back(nd);
		newIdx++;// newIdx is equal to the number of nodes in this module (or network.)
	}

// add outLinks within the newNetwork.
	for (int i = 0; i < numMembers; i++) {
		Node* it = mod->members[i];
		int nid = newNodeID[it->ID()];
		Node* nd_ptr = &(newNetwork.nodes[nid]);
		for (link_iterator link_it = it->outLinks.begin();
				link_it != it->outLinks.end(); link_it++) {
			// check whether the edge within the module or not.
			map<int, int>::iterator m_it = newNodeID.find(link_it->first);
			if (m_it != newNodeID.end()) {
				nd_ptr->outLinks.push_back(
						make_pair(m_it->second, link_it->second));
			}
		}
	}

// add inLinks within the newNetwork based on the generated outLinks above.
	for (vector<Node>::iterator it = newNetwork.nodes.begin();
			it != newNetwork.nodes.end(); it++) {
		for (link_iterator l_it = it->outLinks.begin();
				l_it != it->outLinks.end(); l_it++) {
			newNetwork.nodes[l_it->first].inLinks.push_back(
					make_pair(it->ID(), l_it->second));
		}
	}

	double sum_size_log_size = 0.0;

	for (vector<Node>::iterator it = newNetwork.nodes.begin();
			it != newNetwork.nodes.end(); it++) {
		sum_size_log_size += pLogP(it->Size());
	}
	newNetwork.setAllLogAll(sum_size_log_size);

	for (int i = 0; i < newIdx; i++) {
		newNetwork.modules[i] = Module(i, &newNetwork.nodes[i]);
		newNetwork.nodes[i].setModIdx(i);
	}

	newNetwork.setNModule(newIdx);
	newNetwork.setNNode(newIdx);

	newNetwork.calibrate(1, 1, total_time_calibrate);// This function is run in sequential. also setting tag = 1 for MPISendRecv

	gettimeofday(&endGenNetworkmodule, NULL);

}

void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int> &origNodeID, int numTh, int iteration,
		double& total_time_calibrate) {
	int numMembers = mod->numMembers;
	newNetwork.modules = vector<Module>(numMembers);

	map<int, int> newNodeID;	// key = origNodeID --> value =  newNodeID.

	int newIdx = 0;
	for (vector<Node *>::iterator it = mod->members.begin();
			it != mod->members.end(); it++) {
		newNodeID[(*it)->ID()] = newIdx;
		origNodeID[newIdx] = (*it)->ID();

		Node nd(newIdx, (*it)->Size());
		nd.setNodeWeight((*it)->NodeWeight());
		nd.setTeleportWeight((*it)->TeleportWeight());
		nd.setDanglingSize((*it)->DanglingSize());
		nd.setIsDangling((*it)->IsDangling());

		newNetwork.nodes.push_back(nd);
		newIdx++;// newIdx is equal to the number of nodes in this module (or network.)
	}

	omp_set_num_threads(numTh);	//this statement has nothing to do with my current distributed implementation

	double sum_size_log_size = 0.0;

// add outLinks within the newNetwork.

	for (int i = 0; i < numMembers; i++) {
		Node* it = mod->members[i];
		int nid = newNodeID[it->ID()];
		Node* nd_ptr = &(newNetwork.nodes[nid]);

		for (link_iterator link_it = it->outLinks.begin();
				link_it != it->outLinks.end(); link_it++) {
			// check whether the edge within the module or not.
			map<int, int>::iterator m_it = newNodeID.find(link_it->first);
			if (m_it != newNodeID.end()) {
				nd_ptr->outLinks.push_back(
						make_pair(m_it->second, link_it->second));
			}
		}
	}

// add inLinks within the newNetwork based on the generated outLinks above.
	for (vector<Node>::iterator it = newNetwork.nodes.begin();
			it != newNetwork.nodes.end(); it++) {
		for (link_iterator l_it = it->outLinks.begin();
				l_it != it->outLinks.end(); l_it++) {
			newNetwork.nodes[l_it->first].inLinks.push_back(
					make_pair(it->ID(), l_it->second));
		}
	}

	for (int i = 0; i < numMembers; i++) {
		sum_size_log_size += pLogP(newNetwork.nodes[i].Size());
	}

	newNetwork.setAllLogAll(sum_size_log_size);

	for (int i = 0; i < newIdx; i++) {
		newNetwork.modules[i] = Module(i, &newNetwork.nodes[i]);
		newNetwork.nodes[i].setModIdx(i);
	}

	newNetwork.setNModule(newIdx);
	newNetwork.setNNode(newIdx);

	newNetwork.calibrate(numTh, 2, total_time_calibrate); //also setting tag = 2 for MPISendRecv
}

void print_twoLevel_Cluster(Network network, string networkName,
		string outDir) {
	ofstream outFile;
	ostringstream oss;
	oss << outDir << "/" << networkName << ".tree";
	outFile.open(oss.str().c_str());

	outFile << "# Code length " << network.CodeLength() / log(2.0) << " in "
			<< network.NModule() << " modules." << endl;

	int nModules = network.modules.size();
	int modIdx = 0;
	for (int i = 0; i < nModules; i++) {

		int nMembers = network.modules[i].numMembers;
		if (nMembers > 0)
			modIdx++;

		for (int j = 0; j < nMembers; j++) {
			outFile << modIdx << ":" << j + 1 << " "
					<< network.modules[i].members[j]->Size() << " \""
					<< network.modules[i].members[j]->Name() << "\"" << endl;
		}
	}

	outFile.close();
}

