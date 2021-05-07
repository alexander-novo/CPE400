#include <algorithm>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <queue>
#include <random>
#include <stack>

using json = nlohmann::json;

// See paper for definition of these macros
#define TTL_WINDOW 2.0
#define TTL_BIAS .995

enum Algorithm {
	OSPF,
	TTL,
};

// Command line arguments
struct Arguments {
	json config;
	unsigned seed = 1;
	Algorithm alg;
	bool verbose = false;
};

// A single end of a link between two nodes
struct Connection {
	// nodeIndex is the index of the node on the other end of the link
	// linkConnectionIndex is the index of this connection in the other node's connection list
	// weight is used for routing purposes
	unsigned nodeIndex, linkConnectionIndex;
	double weight, lastTimeWeightUpdated = 0;
	double otherEnergy, otherRate;

	Connection(unsigned nodeIndex, unsigned linkConnectionIndex, double weight)
	    : nodeIndex(nodeIndex), linkConnectionIndex(linkConnectionIndex), weight(weight) {}

	// Operators for the sake of finding max/min
	bool operator>(const Connection& other) const { return weight > other.weight; }
	bool operator<(const Connection& other) const { return weight < other.weight; }
};

struct Node {
	double observationRate, energy, startingEnergy;
	bool edgeNode;
	std::vector<Connection> links;
	std::exponential_distribution<double> dist;
	unsigned shortestLink, observationsMade = 0;
	// A list of observations (timestamp and difference in timestamp since last observation) which have been routed
	// through this node. Always in ascending order of timestamp
	std::vector<std::pair<double, double>> routedObservations;
};

// Network configuration
struct Configuration {
	std::vector<Node> nodes;
	std::vector<unsigned> edges;
	double energyPerTransmit;
};

// Output from the algorithms
struct Output {
	double time = 0;
	unsigned lastObservationNodeIndex, lastRouterNodeIndex, escapedObservations;
};

void simOSPF(const Arguments& arg, Configuration& config, Output& out);
void simTTL(const Arguments& arg, Configuration& config, Output& out);
bool verifyArguments(int argc, char** argv, Arguments& arg, int& err);
void extractConfig(const json& configFile, const Arguments& arg, Configuration& config);
void printHelp();