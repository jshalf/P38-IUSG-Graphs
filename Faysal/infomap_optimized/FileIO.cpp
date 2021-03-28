#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <utility>
#include "Node.h"
#include "Module.h"
#include "FileIO.h"
#include <cstring>
#include <mpi.h>
#include </global/common/sw/cray/cnl7/haswell/metis/5.1.0/intel/19.0.3.199/nbtsmmb/include/metis.h>


using namespace std;

template<class T>
inline string to_string(const T& t) {
	stringstream ss;
	ss << t;
	return ss.str();
}

/*
 * Modified from the Infomap implementation.
 * Reading the given network data in Pajek format.
 */
void load_pajek_format_network(string fName, Network &network) {

	string line;
	string buf;
	string nameBuf;

	/* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
	/* each directed link occurring only once, and link weights > 0.               */
	/* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
	/* Node weights are optional and sets the relative proportion to which         */
	/* each node receives teleporting random walkers. Default value is 1.          */
	/* Example network with three nodes and four directed and weighted links:      */
	/* *Vertices 3                                                                 */
	/* 1 "Name of first node" 1.0                                                  */
	/* 2 "Name of second node" 2.0                                                 */
	/* 3 "Name of third node" 1.0                                                  */
	/* *Arcs 4                                                                     */
	/* 1 2 1.0                                                                     */
	/* 1 3 1.7                                                                     */
	/* 2 3 2.0                                                                     */
	/* 3 2 1.2                                                                     */

	int rank, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cout << "Reading network " << fName << " file\n" << flush;
	ifstream net(fName.c_str());
	network.setNNode(0);

	istringstream ss;
	while (network.NNode() == 0) {
		if (!getline(net, line)) {
			cout << "the network file is not in Pajek format...exiting" << endl;
			exit(-1);
		} else {
			ss.clear();
			ss.str(line);
			ss >> buf;
			if (buf == "*Vertices" || buf == "*vertices"
					|| buf == "*VERTICES") {
				ss >> buf;
				network.setNNode(atoi(buf.c_str()));
			} else {
				cout << "the network file is not in Pajek format...exiting"
						<< endl;
				exit(-1);
			}
		}
	}

	int nNode = network.NNode();

	network.modules = vector<Module>(nNode);
	network.nodes = vector<Node>(nNode);
	network.processTags = new int[8]();
	double totNodeWeights = 0.0;

	// Read node names, assuming order 1, 2, 3, ...
	for (int i = 0; i < nNode; i++) {
		getline(net, line);

		// Assume Vertices format in "id name weight"
		ss.clear();
		ss.str(line);
		ss >> buf;		// read first item, which is ID.
		//network.nodes[i].setID(atoi(buf.c_str()));
		network.nodes[i].setID(i);

		int nameStart = line.find_first_of("\"");
		int nameEnd = line.find_last_of("\"");
		if (nameStart < nameEnd) {
			network.nodes[i].setName(
					string(line.begin() + nameStart + 1,
							line.begin() + nameEnd));
			line = string(line.begin() + nameEnd + 1, line.end());
			ss.clear();
			ss.str(line);
		} else {
			ss.clear();
			ss.str(line);		// After reading ID...
			//ss >> buf;
			//ss >> network.nodeNames[i];
			ss >> buf;
			ss >> buf;
			network.nodes[i].setName(buf);
		}

		buf = "1";
		ss >> buf;
		double nodeWeight = atof(buf.c_str());
		if (nodeWeight <= 0.0)
			nodeWeight = 1.0;
		network.nodes[i].setNodeWeight(nodeWeight);
		totNodeWeights += nodeWeight;
	}

	network.setTotNodeWeights(totNodeWeights);

	// Read the number of links in the network
	getline(net, line);
	ss.clear();
	ss.str(line);
	ss >> buf;

	if (buf != "*Edges" && buf != "*edges" && buf != "*Arcs"
			&& buf != "*arcs") {
		cout << endl << "Number of nodes not matching, exiting" << endl;
		exit(-1);
	}

	double newLinkWeight;
	int NdoubleLinks = 0;

	// Read links in format "from to weight", for example "1 3 0.7"
	while (getline(net, line)) {
		ss.clear();
		ss.str(line);
		ss >> buf;
		int linkEnd1 = atoi(buf.c_str());
		ss >> buf;
		int linkEnd2 = atoi(buf.c_str());
		buf.clear();
		ss >> buf;
		double linkWeight;
		if (buf.empty()) // If no information
			linkWeight = 1.0;
		else
			linkWeight = atof(buf.c_str());
		linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
		linkEnd2--;	// As of now, these values are indices of corresponding nodes in network.nodes vector.

		newLinkWeight = network.Edges[make_pair(linkEnd1, linkEnd2)] +=
				linkWeight;
		if (newLinkWeight > linkWeight)
			NdoubleLinks++;
	}

	net.close();

	network.setNEdge(network.Edges.size());

	cout << "done! (found " << network.NNode() << " nodes and "
			<< network.NEdge() << " edges.)";
	if (NdoubleLinks > 0)
		cout << ", aggregated " << NdoubleLinks
				<< " link(s) defined more than once" << endl;
	else
		cout << endl;
}

void load_csr_format_network(string fName, Network &network, string metisFile) {

	int rank, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char* token;
	int total_vertex = 0, total_edges = 0;
	int start_index = 0;
	int current_vertex = 0;

	char FILENAME[fName.length() + 1];
	strcpy(FILENAME, fName.c_str());

	char METISFILE[metisFile.length() + 1];
	strcpy(METISFILE, metisFile.c_str());

	cout << "Reading network " << fName << " file\n" << flush;
	FILE* fp = fopen(FILENAME, "r");

	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}

	printf("File open successful\n");
	network.setNNode(0);

	int nDoubleLinks = 0;
	double newLinkWeight = 0.0;

	char* line = NULL;
	size_t len = 0;

	if ((getline(&line, &len, fp) != -1)) {
		token = strtok(line, " \n\t");

		if (token != NULL) {
			total_vertex = atoi(token);
			network.setNNode(total_vertex);
			network.modules = vector<Module>(total_vertex);
			network.nodes = vector<Node>(total_vertex);
		}

		token = strtok(NULL, " \n\t");

		if (token != NULL) {
			total_edges = atoi(token);
		}
	}

	for (int i = 0; i < total_vertex; i++) {
		network.nodes[i].setID(i);
		network.nodes[i].setNodeWeight(1.0);
	}

	network.setTotNodeWeights(total_vertex * 1.0);

	while ((getline(&line, &len, fp) != -1)) {
		if (strlen(line) <= 2) {
			printf("end of file reached\n");
			break;
		}

		token = strtok(line, " \n\t");

		while (token != NULL) {
			int edge = atoi(token);
			if (!(network.Edges.count(make_pair(current_vertex, edge)) > 0
					|| network.Edges.count(make_pair(edge, current_vertex)) > 0)) {
				network.Edges[make_pair(current_vertex, edge)] = 1.0;
			}
			token = strtok(NULL, " \n\t");
		}
		current_vertex++;
	}

	network.setNEdge(network.Edges.size());

	network.allNodes = (int*) malloc(total_vertex * sizeof(int));

	if (size == 1) {
		for (int i = 0; i < total_vertex; i++) {
			network.allNodes[i] = 0;
		}
	} else {

		// now read the metis file
		cout << "Reading metis " << metisFile << " file\n" << flush;
		FILE* mp = fopen(METISFILE, "r");

		if (mp == NULL) {
			printf("Could not open metis file\n");
			exit(EXIT_FAILURE);
		}

		while ((getline(&line, &len, mp) != -1)) {
			if (strlen(line) <= 2) {
				printf("end of metis file reached\n");
				break;
			}

			token = strtok(line, " \n\t");
			int index = atoi(token);
			token = strtok(NULL, " \n\t");
			int pRank = atoi(token);
			network.allNodes[index] = pRank;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		fclose(mp);
	}

	fclose(fp);
	if (line) {
		free(line);
	}

	cout << "done! (found " << network.NNode() << " nodes and "
			<< network.NEdge() << " edges.)" << flush;

}

void old_load_csr_format_network(string fName, Network &network) {

	int rank, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	char* token;
	int total_vertex = 0, total_edges = 0;
	int start_index = 0;
	int current_vertex = 0;
	idx_t index = 0;
	int pointer = 0;

	char FILENAME[fName.length() + 1];
	strcpy(FILENAME, fName.c_str());

	cout << "Reading network " << fName << " file\n" << flush;
	FILE* fp = fopen(FILENAME, "r");

	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}

	printf("File open successful\n");
	network.setNNode(0);

	int nDoubleLinks = 0;
	double newLinkWeight = 0.0;

	char* line = NULL;
	size_t len = 0;

	if ((getline(&line, &len, fp) != -1)) {
		token = strtok(line, " \n\t");

		if (token != NULL) {
			total_vertex = atol(token);
			network.setNNode(total_vertex);
			network.modules = vector<Module>(total_vertex);
			network.nodes = vector<Node>(total_vertex);
			network.xadj = (idx_t*) malloc((total_vertex + 1) * sizeof(idx_t));
			network.xadj[index] = 0; //this will be 0 always as this is the starting index
			network.part = (idx_t*) malloc(total_vertex * sizeof(idx_t));
		}

		token = strtok(NULL, " \n\t");

		if (token != NULL) {
			total_edges = atol(token);
			network.adjacency = (idx_t*) malloc(
					(2 * total_edges) * sizeof(idx_t));
		}
	}

	for (int i = 0; i < total_vertex; i++) {
		network.nodes[i].setID(i);
		network.nodes[i].setNodeWeight(1.0);
	}

	network.setTotNodeWeights(total_vertex * 1.0);

	while ((getline(&line, &len, fp) != -1)) {
		if (strlen(line) <= 2) {
			printf("end of file reached\n");
			break;
		}

		token = strtok(line, " \n\t");

		while (token != NULL) {
			int edge = atoi(token);
			if (!(network.Edges.count(make_pair(current_vertex, edge)) > 0
					|| network.Edges.count(make_pair(edge, current_vertex)) > 0)) {
				network.Edges[make_pair(current_vertex, edge)] = 1.0;
			}
			network.adjacency[pointer] = edge;
			pointer++;
			token = strtok(NULL, " \n\t");
		}
		network.vertex_adjacencyIndex.insert(
				make_pair(current_vertex, pointer - 1));
		current_vertex++;
		index++;
		network.xadj[index] = pointer;
	}

	network.setNEdge(network.Edges.size());

	/*	for (int i = 0; i < 2 * network.Edges.size(); i++) {
	 printf("pram rank:%d, network.adjacency[%d]:%d\n", rank, i,
	 network.adjacency[i]);
	 }*/

	network.nvtxs = index;
	network.ncon = 1;
	network.nParts = size;
	/*	int ret = METIS_PartGraphKway(&network.nvtxs, &network.ncon, network.xadj,
	 network.adjacency, NULL, NULL, NULL, &network.nParts, NULL, NULL,
	 NULL, &network.objval, network.part);*/

	fclose(fp);
	if (line) {
		free(line);
	}

	/*	for (unsigned part_i = 0; part_i < network.nvtxs; part_i++) {
	 printf("metis rank:%d, part_i:%d, network.part[%d]:%d\n", rank, part_i,
	 part_i, network.part[part_i]);
	 if (network.part[part_i] == rank) {
	 network.myVertices.push_back(part_i);
	 }
	 }*/
	/*
	 for (int i = 0; i < network.NNode(); i++) {
	 printf("prometis rank:%d, nodeId:%d, processor id:%d\n", rank, i,
	 network.part[i]);
	 }*/

	cout << "done! (found " << network.NNode() << " nodes and "
			<< network.NEdge() << " edges.)" << flush;

}

/*
 *	Read network in the format "FromNodeID	ToNodeID	[Weight]"
 *	not assuming a complete list of nodes, 1..maxnode.
 */
void load_linkList_format_network(string fName, Network &network) {

	string line;
	string buf;
	istringstream ss;

	cout << "Reading network " << fName << " file\n" << flush;
	ifstream net(fName.c_str());
	network.setNNode(0);

	int nDoubleLinks = 0;
	double newLinkWeight = 0.0;

	set<long> Nodes;

	// generate temporary Edge map with long value for large id graphs..
	map<pair<long, long>, double> tempEdges;

	// Read links in format "from to [weight]"
	while (getline(net, line)) {
		if (line[0] != '#') {
			ss.clear();
			ss.str(line);

			ss >> buf;
			long linkEnd1 = atol(buf.c_str());

			ss >> buf;
			long linkEnd2 = atol(buf.c_str());

			buf.clear();
			ss >> buf;
			double linkWeight;
			if (buf.empty()) // If no linkWeight information..
				linkWeight = 1.0;
			else
				linkWeight = atof(buf.c_str());

			// To keep track of all nodes for renaming 0..N-1
			Nodes.insert(linkEnd1);
			Nodes.insert(linkEnd2);

			newLinkWeight = tempEdges[make_pair(linkEnd1, linkEnd2)] +=
					linkWeight;
			if (newLinkWeight > linkWeight)
				nDoubleLinks++;
		}
	}

	net.close();

	network.setNNode(Nodes.size());
	network.setTotNodeWeights(Nodes.size()); // Since no node weight info given, assume uniform weight.
	network.setNEdge(network.Edges.size());

	network.modules = vector<Module>(Nodes.size());
	network.nodes = vector<Node>(Nodes.size());

	// Renaming all nodes if necessary...
	vector<string> nodeNames = vector<string>(Nodes.size());

	long nodeCounter = 0;
	map<long, long> renumber;
	bool renum = false;

	for (set<long>::iterator it = Nodes.begin(); it != Nodes.end(); it++) {
		nodeNames[nodeCounter] = to_string((*it));
		renumber.insert(make_pair((*it), nodeCounter));

		network.nodes[nodeCounter].setID((int) nodeCounter);
		network.nodes[nodeCounter].setName(nodeNames[nodeCounter]);

		if (nodeCounter != (*it))
			renum = true;
		nodeCounter++;
	}

	for (map<pair<long, long>, double>::iterator it = tempEdges.begin();
			it != tempEdges.end(); it++)
		network.Edges.insert(
				make_pair(
						make_pair((int) renumber.find(it->first.first)->second,
								(int) renumber.find(it->first.second)->second),
						it->second));

	cout << "done! (found " << network.NNode() << " nodes and "
			<< network.NEdge() << " edges.)";
	if (nDoubleLinks > 0)
		cout << ", aggregated " << nDoubleLinks
				<< " link(s) defined more than once" << endl;
	else
		cout << endl;

	cout << "done! (found " << network.nodes.size() << " nodes and "
			<< network.Edges.size() << " edges.)" << endl;
}

void write_cluster_result(string outFile, Network &finalNetwork) {

}
