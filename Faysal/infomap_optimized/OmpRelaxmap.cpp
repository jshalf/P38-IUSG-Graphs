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
#include "global.h"
#include </opt/intel/vtune_profiler_2020.2.0.610396/include/ittnotify.h>
#include <chrono>


using namespace std;

double totalExecutionTime;
double convert2SpNodeTime;
double computePageRankTime;
double prioritizeMoveTime;
double prioritizeSpNodeMoveTime;
double updateMembersTime;
double initiationTime;
double networkReadTime;
double calibrationTime;
double prioritizeMoveNodeAccessTime;
double prioritizeMoveNodeAccessTimeChunk;
double prioritizeMoveModuleAccessTime;
double prioritizeSpMoveNodeAccessTime;
double prioritizeSpMoveModuleAccessTime;
double randomNumberGeneratorTime;
double bestModuleSelectionTime;
double prioritizeSpNodeModuleAccessTimeChunk;
double bestSpModuleSelectionTime;
vector<double>threadtimes(68, 0.0);
double prioritizeSpNodeCommunicationTime;
double prioritizeNodeCommunicationTime;
double prioritizeSpNodeSyncTime;
double prioritizeNodeSyncTime;
double convertModulePart1;
double convertModulePart2;
double convertModulePart3;

int iterations_updateMembers;
int iteration_convertModules;
int iteration_prioritizeMove;
int iteration_prioritizeSpMove;


unsigned stou(char *s) {
	return strtoul(s, (char **) NULL, 10);
}

void stochastic_greedy_partition(Network &network, int numTh, double threshold,
		double vThresh, int maxIter, bool prior, bool fineTune, bool fast,
		bool inLoop);
void partition_module_network(Network &network, int numTh, double threshold,
		int maxIter, bool fast);
void generate_sub_modules(Network &network, int numTh, double threshold,
		int maxIter);
void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int>& origNodeID, int numTh, int iteration);
void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int> &origNodeID, int iteration);

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
	
	//__itt_pause();
	totalExecutionTime = 0.0;
	convert2SpNodeTime = 0.0;
	computePageRankTime = 0.0;
	prioritizeMoveTime = 0.0;
	prioritizeSpNodeMoveTime = 0.0;
	updateMembersTime = 0.0;
	networkReadTime = 0.0;
	initiationTime = 0.0;
	calibrationTime = 0.0;
	prioritizeMoveNodeAccessTime = 0.0;
	prioritizeMoveModuleAccessTime = 0.0;
	prioritizeSpMoveNodeAccessTime = 0.0;
	prioritizeSpMoveModuleAccessTime = 0.0;
	prioritizeMoveNodeAccessTimeChunk = 0.0;
	randomNumberGeneratorTime = 0.0;
	bestModuleSelectionTime = 0.0;
	prioritizeSpNodeModuleAccessTimeChunk = 0.0;
	bestSpModuleSelectionTime = 0.0;
	prioritizeSpNodeCommunicationTime = 0.0;
	prioritizeNodeCommunicationTime = 0.0;
	prioritizeSpNodeSyncTime = 0.0;
	prioritizeNodeSyncTime = 0.0;
	convertModulePart1 = 0.0;
	convertModulePart2 = 0.0;
	convertModulePart3 = 0.0;


	iterations_updateMembers = 0;
	iteration_convertModules = 0;
	iteration_prioritizeMove = 0;
	iteration_prioritizeSpMove = 0;



	string outDir = string(argv[8]);
	int maxIter = atoi(argv[7]);	// Set the maximum number of iteration..
	int Ntrials = atoi(argv[4]);  // Set number of partition attempts
	int numThreads = atoi(argv[3]);	// Set number of threads...
	string line;
	string buf;

	MTRand *R = new MTRand(stou(argv[1]));

	string infile = string(argv[2]);
	string networkFile = string(argv[2]);
	string networkName(networkFile.begin() + networkFile.find_last_of("/"), networkFile.begin() + networkFile.find_last_of("."));
	string networkType(infile.begin() + infile.find_last_of("."), infile.end());

	double threshold = atof(argv[5]);
	double vThresh = atof(argv[6]);	// vertex-threshold: threshold for each vertex-movement.

	//cout << "Threshold = " << threshold << ", Vertex-Threshold = " << vThresh << endl;

	vThresh *= -1;	// change the threshold value for negative.

	string priorFlag = string(argv[9]);

	bool prior = false;
	if (priorFlag == "prior")
		prior = true;

	string metisFile = string(argv[10]);

	bool includeSelfLinks = false;
	if (argc == 12) 
	{
		string selfLinks(argv[11]);
		if (selfLinks == "selflinks")
		{
			includeSelfLinks = true;
		}
	}

	Network origNetwork;	// default constructor is called.

	origNetwork.R = R;

	// time values for measuring elapsed times for each step...
	struct timeval allStart, allEnd;
	struct timeval noIOstart, noIOend;
	struct timeval start, end;

	gettimeofday(&allStart, NULL);
	gettimeofday(&start, NULL);

    	auto tik = std::chrono::high_resolution_clock::now();

	if (networkType == ".net") 
	{
		load_pajek_format_network(networkFile, origNetwork);
	} 
	else if (networkType == ".csr") 
	{
		load_csr_format_network(networkFile, origNetwork, metisFile);
	} 
	else 
	{
		load_linkList_format_network(networkFile, origNetwork);
	}
		
	
	gettimeofday(&end, NULL);
	
	cout << "Time for reading input data : " << elapsedTimeInSec(start, end) << " (sec)" << endl;

	auto ding = std::chrono::high_resolution_clock::now();

	networkReadTime = std::chrono::duration_cast<std::chrono::nanoseconds>(ding - tik).count();

	gettimeofday(&noIOstart, NULL);

	int nNode = origNetwork.NNode();

	double totNodeWeights = origNetwork.TotNodeWeights();

	cout << "total Node weights = " << totNodeWeights << endl;

	//here let's build the parmetis graph
	//origNetwork.buildParMetis();

	gettimeofday(&start, NULL);

	for (int i = 0; i < nNode; i++) 
	{
		origNetwork.nodes[i].setTeleportWeight(origNetwork.nodes[i].NodeWeight() / totNodeWeights);
	}

	int NselfLinks = 0;

	for (map<pair<int, int>, double>::iterator it = origNetwork.Edges.begin(); it != origNetwork.Edges.end(); it++) 
	{

		int from = it->first.first;
		int to = it->first.second;
		double weight = it->second;
		if (weight > 0.0) 
		{
			if (from == to) 
			{
				NselfLinks++;
				//if(includeSelfLinks)
				//	origNetwork.nodes[from]->selfLink += weight;
			} 
			else 
			{
				origNetwork.nodes[from].outLinks.push_back(make_pair(to, weight));
				// we will going to update inLinks, after we got final flow of the network.
				//origNetwork.nodes[to].inLinks.push_back(make_pair(from,weight));
			}
		}
	}

	if (includeSelfLinks)
	//cout << "including " <<  NselfLinks << " self link(s)." << endl;
	{
		cout << "current version always excludes self links.\nignoring " << NselfLinks << " self link(s)." << endl;
	} 
	else 
	{
		cout << "ignoring " << NselfLinks << " self link(s)." << endl;
	}

	//Swap vector to free memory
	map<pair<int, int>, double>().swap(origNetwork.Edges);

	cout << "DONE: Parsing the given network  ..." << endl;

	gettimeofday(&end, NULL);

	cout << "Time for parsing the given network : " << elapsedTimeInSec(start, end) << " (sec)" << endl;

	gettimeofday(&start, NULL);
	// Initialization.

	origNetwork.initiate(numThreads);

	// Now update inLinks..
	for (int i = 0; i < nNode; i++) 
	{
		int nOutLinks = origNetwork.nodes[i].outLinks.size();

		for (int j = 0; j < nOutLinks; j++)
		{
			origNetwork.nodes[origNetwork.nodes[i].outLinks[j].first].inLinks.push_back(make_pair(i, origNetwork.nodes[i].outLinks[j].second));
		}

	}

	gettimeofday(&end, NULL);

	printf("DONE: Initiate() ... for rank:%d in %f sec\n", rank, elapsedTimeInSec(start, end));
	cout << "Initial Code Length: " << origNetwork.CodeLength() / log(2.0) << " in " << origNetwork.NModule() << " modules.\n";

	// copy size of each node for print in order.
	vector<double> nodeSize(nNode);

	for (int i = 0; i < nNode; i++) 
	{
		nodeSize[i] = origNetwork.nodes[i].Size();
	}

	cout << "Now partition the network starts...\n";

	gettimeofday(&start, NULL);

	bool fineTune = true;
	bool fast = false;	// This will be true only for sub-module partitioning...

	int step = 1;

	// Initial SuperStep running ...
	double oldCodeLength = origNetwork.CodeLength();

	stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh, maxIter, prior, fineTune, fast, false);

	printf("SuperStep[%d] - codeLength = %f in %d modules.\n", step, origNetwork.CodeLength() / log(2.0), origNetwork.NModule());

	bool nextIter = true;

	if ((oldCodeLength - origNetwork.CodeLength()) / log(2.0) < threshold || origNetwork.CodeLength() <= 0.0) 
	{
		nextIter = false;
	}

	struct timeval subStart, subEnd;

	//nextIter = false; // I am forcefully doing it to stop executing the while loop below for now

	while (nextIter) {

		oldCodeLength = origNetwork.CodeLength();

		stochastic_greedy_partition(origNetwork, numThreads, threshold, vThresh, maxIter, prior, fineTune, fast, true);

		step++;

		printf("SuperStep[%d] - codeLength = %f in %d modules.\n", step, origNetwork.CodeLength() / log(2.0), origNetwork.NModule());

		if ((oldCodeLength - origNetwork.CodeLength()) / log(2.0) < threshold || origNetwork.CodeLength() < 0.0) 
		{
			nextIter = false;			
		}

		//TODO: for now I am ignoring the operation of reconstructing network based on found modules

		/*fineTune = !fineTune; // fine-tune and coarse-tune will be done alternatively.

		 if (nextIter && !fineTune) {
		 // Next iteration will be Coarse Tune.
		 gettimeofday(&subStart, NULL);

		 generate_sub_modules(origNetwork, numThreads, threshold, maxIter);

		 gettimeofday(&subEnd, NULL);
		 printf("Time for finding sub-modules for rank:%d is:%f sec\n", rank,
		 elapsedTimeInSec(subStart, subEnd));
		 }*/
	}

	gettimeofday(&end, NULL);
	cout << "Time for partitioning : " << elapsedTimeInSec(start, end) << " (sec)" << endl;

	cout << "DONE: Code Length = " << origNetwork.CodeLength() / log(2.0) << " in ";
	cout << origNetwork.NModule() << " modules, with " << nNode << " nodes.\n" << endl;

	gettimeofday(&noIOend, NULL);
	gettimeofday(&allEnd, NULL);

	cout << "Overall Elapsed Time for Module Detection (w/o file IO): " << elapsedTimeInSec(noIOstart, noIOend) << " (sec)" << endl;
	cout << "Overall Elapsed Time for Module Detection (w/ file Reading): " << elapsedTimeInSec(allStart, allEnd) << " (sec)" << endl;

	cout << "\nComputed Code Length = "<< origNetwork.calculateCodeLength() / log(2.0) << endl;

	double modularity = origNetwork.calculateModularityScore();

	printf("\nmodularity score for the network:%f\n", modularity);

	double conductance = origNetwork.calculateConductance();

	printf("\nconductance of the network:%f\n\n", conductance);


//Print two-level clustering result in .tree file
	print_twoLevel_Cluster(origNetwork, networkName, outDir);

// Print partition in Pajek's .clu format
	ofstream outFile;
	ostringstream oss;

	oss.str("");
	oss << outDir << "/" << networkName << ".clu";
	outFile.open(oss.str().c_str());
	outFile << "*Vertices " << nNode << "\x0D\x0A";
	for (int i = 0; i < nNode; i++) 
	{
		outFile << origNetwork.nodes[i].ModIdx() + 1 << "\x0D\x0A";
	}
	outFile.close();
  
	auto tok = std::chrono::high_resolution_clock::now();
	totalExecutionTime = std::chrono::duration_cast<std::chrono::nanoseconds>(tok - tik).count();

     	int myrank;

    	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    	if(myrank==0)
	{
		printf("========totalexecutiontime:%f========\n", totalExecutionTime*(1e-9));
		printf("========networkReadTime:%f========\n", networkReadTime*(1e-9));
		printf("========initiationTime:%f========\n", initiationTime*(1e-9));
		printf("========calibrationTime:%0.9f========\n", calibrationTime*(1e-9));
		printf("========computePageRankTime:%f========\n", computePageRankTime*(1e-9));
		printf("========convert2SpNodeTime:%0.9f========in iterations:%d========\n", convert2SpNodeTime*(1e-9), iteration_convertModules);
		printf("========convertModulePart1:%0.9f========\n", convertModulePart1*(1e-9));
		printf("========convertModulePart2:%0.9f========\n", convertModulePart2*(1e-9));
		printf("========convertModulePart3:%0.9f========\n", convertModulePart3*(1e-9));
		printf("========updateMembersTime:%0.9f========in %d iterations========\n", updateMembersTime*(1e-9), iterations_updateMembers);
		printf("========prioritizeMoveTime:%0.9f========in %d iterations========\n", prioritizeMoveTime*(1e-9), iteration_prioritizeMove);
		printf("========prioritizeSpNodeMoveTime:%0.9f========in %d iterations========\n", prioritizeSpNodeMoveTime*(1e-9), iteration_prioritizeSpMove);
		printf("========prioritizeMoveNodeAccessTimeChunk:%0.9f=======bestModuleSelectionTime:%0.9f==========\n", prioritizeMoveNodeAccessTimeChunk*(1e-9), bestModuleSelectionTime*(1e-9));
		printf("========prioritizeSpNodeModuleAccessTimeChunk:%0.9f=======bestSpModuleSelectionTime:%0.9f==========\n", prioritizeSpNodeModuleAccessTimeChunk*(1e-9), bestSpModuleSelectionTime*(1e-9));
		printf("========prioritizeNodeCommunicationTime:%0.9f========\n", prioritizeNodeCommunicationTime*(1e-9));
		printf("========prioritizeSpNodeCommunicationTime:%0.9f========\n", prioritizeSpNodeCommunicationTime*(1e-9));
		printf("========prioritizeNodeSyncTime:%0.9f========\n", prioritizeNodeSyncTime*(1e-9));
		printf("========prioritizeSpNodeSyncTime:%0.9f========\n", prioritizeSpNodeSyncTime*(1e-9));
		
		/*
		for(int k = 0; k < 68; k++)
		{
			printf("threadtimes[%d]:%0.9f\n", k, threadtimes[k]*(1e-9));
		}
		*/		
	}
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
void stochastic_greedy_partition(Network &network, int numTh, double threshold, double vThresh, int maxIter, bool prior, bool fineTune, bool fast, bool inLoop) 
{

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

	for (int i = 0; i < nActiveUnits; i++) 
	{
		network.activeNodes[i] = i;
		network.isActives[i] = 0;	// initially set inactive nodes.
		network.activeVertices.insert(make_pair(i, 1));
	}

	int numMoved = 0;

	while (!stop && iter < maxIter) 
	{
		gettimeofday(&inner_T1, NULL);

		oldCodeLength = network.CodeLength();

		if (fineTune) 
		{
			if (prior) 
			{
				numMoved = network.prioritize_move(vThresh, numTh, inLoop);
			} 
			else 
			{
				numMoved = network.move(numTh);
			}
		} 
		else 
		{
			if (prior) 
			{
				numMoved = network.prioritize_moveSPnodes(vThresh, numTh, iter, inLoop);
			} 
			else 
			{
				numMoved = network.moveSuperNodes(numTh);// If at least one node is moved, return true. Otherwise, return false.
			}
		}
		iter++;

		if (network.CodeLength() <= 0 || ((oldCodeLength - network.CodeLength()) / log(2.0) < threshold)) 
		{
			stop = true;	//moved = false;
		}
		gettimeofday(&inner_T2, NULL);
	}

	int outerLoop = 1;

	network.updateMembersInModule();

	tSequential += elapsedTimeInSec(seq_T1, seq_T2);

	if (fast)
		return;

	double tConvert = 0.0;

	do 
	{
		oldCodeLength = network.CodeLength();
		stop = false;

		network.convertModulesToSuperNodes(numTh);

		nActiveUnits = network.superNodes.size();

		// set initial active nodes list ...
		vector<char>(nActiveUnits).swap(network.isActives);
		vector<int>(nActiveUnits).swap(network.activeNodes);
		for (int i = 0; i < nActiveUnits; i++) {
			network.activeNodes[i] = i;	// initially all active units are active.
			network.isActives[i] = 0;	// initially set inactive nodes.
		}

		int spIter = 0;
		while (!stop && spIter < maxIter) 
		{

			double innerOldCodeLength = network.CodeLength();

			if (prior) 
			{
				numMoved = network.prioritize_moveSPnodes(vThresh, numTh, spIter, inLoop);
			} 
			else
			{
				numMoved = network.moveSuperNodes(spIter);
			}
			spIter++;

			
			if (network.CodeLength() <= 0.0 || (innerOldCodeLength - network.CodeLength()) / log(2.0) < threshold) 
			{
				stop = true;
			}
			gettimeofday(&inner_T2, NULL);
		}

		network.updateMembersInModule();

		outerLoop++;

		tag++;

	} while ((oldCodeLength - network.CodeLength()) / log(2.0) > threshold);

}

/**
 *	This function will be called for partitioning sub-Module of each module of the original graph.
 *	Thus, we would like to reduce printing from this function for providing high-level log.
 */
void partition_module_network(Network &network, int numTh, double threshold, int maxIter, bool fast) {

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

	while (!stop && iter < maxIter) 
	{
		oldCodeLength = network.CodeLength();

		if (numTh == 1) {
			numMoved = network.move(iter);
		} else {
			numMoved = network.parallelMove(numTh, tSequential);
		}

		iter++;

		if ((oldCodeLength - network.CodeLength()) / log(2.0) < threshold)
			stop = true;
	}

	int outerLoop = 1;

	network.updateMembersInModule();

	if (fast)
	{
		return;
	}

	do {
		oldCodeLength = network.CodeLength();

		stop = false;
		network.convertModulesToSuperNodes(numTh);

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

		network.updateMembersInModule();

		outerLoop++;

	} while ((oldCodeLength - network.CodeLength()) / log(2.0) > threshold);

	gettimeofday(&endPartition, NULL);

	//printf("time for partition_module_network in rank:%d is %f\n", rank,

}

void generate_sub_modules(Network &network, int numTh, double threshold, int maxIter) {

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

			generate_network_from_module(newNetwork, mod, origNodeID, iteration);
			newNetwork.R = Rand;

			partition_module_network(newNetwork, 1, threshold, maxIter, true);	// fast = true..

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
		generate_network_from_module(newNetwork, mod, origNodeID, numTh, iteration);
		newNetwork.R = Rand;

		partition_module_network(newNetwork, numTh, threshold, maxIter, true); // fast = true..

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
		map<int, int> &origNodeID, int iteration) 
{

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

	newNetwork.calibrate(1, 1);// This function is run in sequential. also setting tag = 1 for MPISendRecv

	gettimeofday(&endGenNetworkmodule, NULL);

}

void generate_network_from_module(Network &newNetwork, Module* mod,
		map<int, int> &origNodeID, int numTh, int iteration) {
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

	newNetwork.calibrate(numTh, 2); //also setting tag = 2 for MPISendRecv
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

