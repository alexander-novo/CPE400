#include "main.h"

int main(int argc, char** argv) {
	int err;
	Arguments arg;

	// Read in command line arguments
	if (!verifyArguments(argc, argv, arg, err)) { return err; }

	// Read in network configuration
	Configuration config;
	try {
		extractConfig(arg.config, arg, config);
	} catch (json::exception e) {
		std::cout << "Given configuration file was invalid:\n" << e.what() << std::endl;
		return 2;
	}

	// Complete sim depending on algorithm chosen
	Output out;
	switch (arg.alg) {
		case Algorithm::OSPF:
			simOSPF(arg, config, out);
			break;
		case Algorithm::TTL:
			simTTL(arg, config, out);
			break;
	}

	// Output results
	std::cout << "Collected " << out.escapedObservations << " observations in " << out.time << "s. Stopped after node "
	          << out.lastObservationNodeIndex << " made an observation and attempted to route through node "
	          << out.lastRouterNodeIndex << ", which did not have enough energy to route it.\n";

	std::cout << "Node   Remaining energy   Observations made   Observations routed\n";
	for (unsigned i = 0; i < config.nodes.size(); i++) {
		const Node& node = config.nodes[i];

		printf("%4u   %16f   %17d   %19u\n", i, node.energy, node.observationsMade,
		       (unsigned) round(std::max(
		           0., (node.startingEnergy - node.energy) / config.energyPerTransmit - node.observationsMade)));
	}
}

void simOSPF(const Arguments& arg, Configuration& config, Output& out) {
	if (arg.verbose) std::cout << "Edges: ";
	// Start by filling out routing tables - beginning from each edge / escape node, since they have a minimum weight
	// link (going to the external gateway)
	for (unsigned edgeIndex : config.edges) {
		std::priority_queue<std::tuple<unsigned, unsigned, unsigned>> queue;
		Node& edge = config.nodes[edgeIndex];
		for (Connection con : edge.links) { queue.emplace(1, con.nodeIndex, con.linkConnectionIndex); }

		if (arg.verbose) std::cout << edgeIndex << ' ';

		while (!queue.empty()) {
			auto [weight, nodeIndex, connectionIndex] = queue.top();
			Node& node                                = config.nodes[nodeIndex];

			queue.pop();

			if (node.links[connectionIndex].weight <= weight) { continue; }

			node.links[connectionIndex].weight = weight;

			for (Connection con : node.links) {
				if (con.nodeIndex != nodeIndex) queue.emplace(weight + 1, con.nodeIndex, con.linkConnectionIndex);
			}
		}
	}
	if (arg.verbose) std::cout << '\n';

	// Calculate the shortest path - the connection with the lowest weight - from each node
	for (Node& node : config.nodes) {
		node.shortestLink = std::distance(
		    node.links.begin(), std::max_element(node.links.begin(), node.links.end(), std::greater<Connection>()));
	}

	// Print out connection status if verbose mode was turned on
	if (arg.verbose) {
		std::cout << "Connections:\n  Nodes      Weights\n";
		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			for (unsigned j = 0; j < node.links.size(); j++) {
				const Connection& con = node.links[j];
				if (con.nodeIndex > i) {
					std::cout << std::setw(2) << i << " <-> " << std::setw(2) << con.nodeIndex << "    " << con.weight
					          << " <-> " << config.nodes[con.nodeIndex].links[con.linkConnectionIndex].weight << '\n';
				}
			}
		}

		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			std::cout << "Node " << std::setw(2) << i << " shortest path: " << std::setw(2)
			          << node.links[node.shortestLink].nodeIndex << '\n';
		}
	}

	std::vector<std::pair<double, unsigned>> observationEvents;
	std::mt19937 gen(arg.seed);

	// Setup initial observations
	for (unsigned i = 0; i < config.nodes.size(); i++) {
		Node& node = config.nodes[i];
		observationEvents.emplace_back(node.dist(gen), i);
	}

	// Make min heap so that the next observation is always in the front
	std::make_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

	if (arg.verbose) {
		std::cout << "Initial Observations: (heap)\n"
		          << "Time       Node\n";
		for (const std::pair<double, unsigned>& observation : observationEvents) {
			printf("%8f   %4d\n", observation.first, observation.second);
		}
	}

	// Start simulation
	out.time                = 0;
	out.escapedObservations = 0;
	while (true) {
		// Get next observation event, maintain heap
		std::pair<double, unsigned> observation = observationEvents.front();
		Node& observingNode                     = config.nodes[observation.second];
		std::pop_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

		// Step through time
		double deltaT = observation.first;
		out.time += deltaT;
		std::for_each(observationEvents.begin(), observationEvents.end(),
		              [deltaT](std::pair<double, unsigned>& observation) { observation.first -= deltaT; });

		// Add next observation event for this node
		observationEvents.back() = {observingNode.dist(gen), observation.second};
		std::push_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());
		observingNode.observationsMade++;

		// Route
		unsigned currentNodeIndex = observation.second;
		while (!config.nodes[currentNodeIndex].edgeNode &&
		       config.nodes[currentNodeIndex].energy >= config.energyPerTransmit) {
			Node& currentNode = config.nodes[currentNodeIndex];

			currentNode.energy -= config.energyPerTransmit;
			currentNodeIndex = currentNode.links[currentNode.shortestLink].nodeIndex;
		}

		// If we were able to route to an edge node, we collected the data.
		// Otherwise, stop simulation
		Node& lastNode = config.nodes[currentNodeIndex];
		if (lastNode.edgeNode && lastNode.energy >= config.energyPerTransmit) {
			lastNode.energy -= config.energyPerTransmit;
			out.escapedObservations++;
		} else {
			out.lastObservationNodeIndex = observation.second;
			out.lastRouterNodeIndex      = currentNodeIndex;
			break;
		}
	}
}

void simTTL(const Arguments& arg, Configuration& config, Output& out) {
	if (arg.verbose) std::cout << "Edges: ";
	// Start by filling out routing tables - beginning from each edge / escape node, since their TTL will always depend
	// only on their own parameters (they have an infinite TTL neighbor - the gateway).
	for (unsigned edgeIndex : config.edges) {
		std::priority_queue<std::tuple<double, double, unsigned, unsigned>> queue;
		Node& edge = config.nodes[edgeIndex];
		for (Connection& con : edge.links) {
			queue.emplace(edge.energy, edge.observationRate, con.nodeIndex, con.linkConnectionIndex);
		}

		if (arg.verbose) std::cout << edgeIndex << ' ';

		while (!queue.empty()) {
			auto [energy, observationRate, nodeIndex, connectionIndex] = queue.top();
			Node& node                                                 = config.nodes[nodeIndex];

			queue.pop();

			// TTL and TTL' (see paper)
			double ttl1 = energy / config.energyPerTransmit / (observationRate + node.observationRate);
			double ttl2 = node.energy / config.energyPerTransmit / node.observationRate;
			double ttl  = std::min(ttl1, ttl2);

			Connection& con = node.links[connectionIndex];
			if (con.weight >= ttl) { continue; }

			// otherEnergy = e' and otherRate = lambda' (see paper)
			con.weight      = ttl;
			con.otherEnergy = energy;
			con.otherRate   = observationRate;

			for (Connection con : node.links) {
				if (con.nodeIndex != nodeIndex) {
					if (ttl1 < ttl2) {
						queue.emplace(energy, observationRate + node.observationRate, con.nodeIndex,
						              con.linkConnectionIndex);
					} else {
						queue.emplace(node.energy, node.observationRate, con.nodeIndex, con.linkConnectionIndex);
					}
				}
			}
		}
	}
	if (arg.verbose) std::cout << '\n';

	// Fill out routing tables, picking the neighbor with the highest TTL
	for (Node& node : config.nodes) {
		node.shortestLink = std::distance(node.links.begin(), std::max_element(node.links.begin(), node.links.end()));
	}

	// Print connection information if verbose. Note - routing table should be identical to OSPF at this point
	if (arg.verbose) {
		std::cout << "Connections:\n  Nodes      Weights\n";
		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			for (unsigned j = 0; j < node.links.size(); j++) {
				const Connection& con = node.links[j];
				if (con.nodeIndex > i) {
					std::cout << std::setw(2) << i << " <-> " << std::setw(2) << con.nodeIndex << "    " << con.weight
					          << " <-> " << config.nodes[con.nodeIndex].links[con.linkConnectionIndex].weight << '\n';
				}
			}
		}

		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			std::cout << "Node " << std::setw(2) << i << " shortest path: " << std::setw(2)
			          << node.links[node.shortestLink].nodeIndex << '\n';
		}
	}

	// generate initial observation events
	std::vector<std::pair<double, unsigned>> observationEvents;
	std::mt19937 gen(arg.seed);

	for (unsigned i = 0; i < config.nodes.size(); i++) {
		Node& node = config.nodes[i];
		observationEvents.emplace_back(node.dist(gen), i);
	}

	// Make min heap
	std::make_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

	if (arg.verbose) {
		std::cout << "Initial Observations: (heap)\n"
		          << "Time       Node\n";
		for (const std::pair<double, unsigned>& observation : observationEvents) {
			printf("%8f   %4d\n", observation.first, observation.second);
		}
	}

	// Start simulation
	out.time                = 0;
	out.escapedObservations = 0;
	while (true) {
		// Get next observation event, maintain heap
		std::pair<double, unsigned> observation = observationEvents.front();
		Node& observingNode                     = config.nodes[observation.second];
		std::pop_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

		// Step through time - this is big delta t (see paper)
		double deltaT = observation.first;
		out.time += deltaT;
		std::for_each(observationEvents.begin(), observationEvents.end(),
		              [deltaT](std::pair<double, unsigned>& observation) { observation.first -= deltaT; });

		// Generate new observation event from this node
		observationEvents.back() = {observingNode.dist(gen), observation.second};
		std::push_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());
		observingNode.observationsMade++;

		// Route - keep track of the path we routed through to update TTL
		unsigned currentNodeIndex = observation.second;
		std::stack<unsigned> routedPath;
		while (!config.nodes[currentNodeIndex].edgeNode &&
		       config.nodes[currentNodeIndex].energy >= config.energyPerTransmit) {
			Node& currentNode = config.nodes[currentNodeIndex];

			currentNode.energy -= config.energyPerTransmit;
			routedPath.push(currentNodeIndex);
			currentNodeIndex = currentNode.links[currentNode.shortestLink].nodeIndex;

			// This is little delta t (see paper)
			double deltaTSinceLastRoute =
			    config.nodes[currentNodeIndex].routedObservations.size() > 0 ?
                    out.time - config.nodes[currentNodeIndex].routedObservations.back().first :
                    out.time;

			config.nodes[currentNodeIndex].routedObservations.emplace_back(out.time, deltaTSinceLastRoute);
		}

		// If we made it to an edge node - we collected an observation and need to update TTL
		Node& lastNode = config.nodes[currentNodeIndex];
		if (lastNode.edgeNode && lastNode.energy >= config.energyPerTransmit) {
			lastNode.energy -= config.energyPerTransmit;
			out.escapedObservations++;

			double lastEnergy;
			double lastRate;

			// Go back through path and update TTL in reverse order
			while (!routedPath.empty()) {
				Node& currentNode = config.nodes[currentNodeIndex];
				currentNodeIndex  = routedPath.top();
				routedPath.pop();
				Node& prevNode = config.nodes[currentNodeIndex];

				// This is updating lambda based on the window
				double perceivedRouteReceiveRate =
				    std::distance(
				        currentNode.routedObservations.crbegin(),
				        std::find_if(
				            currentNode.routedObservations.crbegin(), currentNode.routedObservations.crend(),
				            [&out](const std::pair<double, int>& obs) { return out.time - obs.first < TTL_WINDOW; })) /
				    std::min(out.time, TTL_WINDOW);

				// TTL and TTL' (see paper)
				double ttl1 = currentNode.energy / config.energyPerTransmit /
				              (perceivedRouteReceiveRate + currentNode.observationRate);
				double ttl2 = currentNode.edgeNode ?
                                  std::numeric_limits<double>::infinity() :
                                  lastEnergy / config.energyPerTransmit /
				                      (perceivedRouteReceiveRate + currentNode.observationRate + lastRate);

				// Small delta t again
				double deltaTSinceLastRoute = currentNode.routedObservations.back().second;

				// Updating old TTL information from neighbors, using bias (see paper)
				std::for_each(
				    prevNode.links.begin(), prevNode.links.end(), [deltaTSinceLastRoute, &out](Connection& con) {
					    con.weight -= std::min(
					        con.weight, pow(TTL_BIAS, out.time - con.lastTimeWeightUpdated) * deltaTSinceLastRoute);
				    });

				// Updating routing table
				prevNode.links[prevNode.shortestLink].weight                = std::min(ttl1, ttl2);
				prevNode.links[prevNode.shortestLink].lastTimeWeightUpdated = out.time;
				prevNode.shortestLink                                       = std::distance(prevNode.links.begin(),
                                                      std::max_element(prevNode.links.begin(), prevNode.links.end()));

				// Since we may have changed the shortest link, we need to recalc ttl using what we think we know about
				// the shortest link
				ttl2 = currentNode.edgeNode ?
                           std::numeric_limits<double>::infinity() :
                           prevNode.links[prevNode.shortestLink].otherEnergy / config.energyPerTransmit /
				               (perceivedRouteReceiveRate + currentNode.observationRate +
				                prevNode.links[prevNode.shortestLink].otherRate);

				// Propagation information e' and lambda'
				if (ttl1 < ttl2) {
					lastEnergy = currentNode.energy;
					lastRate   = perceivedRouteReceiveRate + currentNode.observationRate;
				} else {
					lastEnergy = prevNode.links[prevNode.shortestLink].otherEnergy;
					lastRate   = perceivedRouteReceiveRate + currentNode.observationRate +
					           prevNode.links[prevNode.shortestLink].otherRate;
				}
			}
		} else {
			out.lastObservationNodeIndex = observation.second;
			out.lastRouterNodeIndex      = currentNodeIndex;
			break;
		}
	}

	// Final state information for debugging.
	if (arg.verbose) {
		std::cout << "Final Weights:\n  Nodes      Weights\n";
		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			for (unsigned j = 0; j < node.links.size(); j++) {
				const Connection& con = node.links[j];
				if (con.nodeIndex > i) {
					std::cout << std::setw(2) << i << " <-> " << std::setw(2) << con.nodeIndex << "    " << con.weight
					          << " <-> " << config.nodes[con.nodeIndex].links[con.linkConnectionIndex].weight << '\n';
				}
			}
		}

		for (unsigned i = 0; i < config.nodes.size(); i++) {
			const Node& node = config.nodes[i];
			std::cout << "Node " << std::setw(2) << i << " shortest path: " << std::setw(2)
			          << node.links[node.shortestLink].nodeIndex << '\n';
		}
	}
}

bool verifyArguments(int argc, char** argv, Arguments& arg, int& err) {
	if (argc < 2 || (argc < 3 && strcmp(argv[1], "-h") && strcmp(argv[1], "--help"))) {
		std::cout << "Missing operand.\n\n";
		err = 1;
		printHelp();
		return false;
	}

	if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
		printHelp();
		return false;
	}

	// Required arguments
	if (!strcmp(argv[1], "ospf")) {
		arg.alg = Algorithm::OSPF;
	} else if (!strcmp(argv[1], "ttl")) {
		arg.alg = Algorithm::TTL;
	} else {
		std::cout << "Algorithm type \"" << argv[1] << "\" not recognised.\n\n";
		err = 1;
		printHelp();
		return false;
	}
	std::ifstream configFile(argv[2]);
	if (!configFile.is_open() || !configFile) {
		std::cout << "Could not open file \"" << argv[2] << "\".\n";
		err = 1;
		return false;
	}
	try {
		// Ignore comments so that we can comment the config files
		arg.config = json::parse(configFile, nullptr, true, true);
	} catch (json::exception e) { std::cout << "Could not parse \"" << argv[1] << "\" as JSON.\n"; }

	for (unsigned i = 3; i < argc; i++) {
		if (!strcmp(argv[i], "-s")) {
			if (i + 1 >= argc) {
				std::cout << "Missing seed argument.\n\n";
				err = 1;
				printHelp();
				return false;
			}

			char* endptr;
			arg.seed = strtol(argv[i + 1], &endptr, 10);
			if (endptr == argv[i + 1]) {
				std::cout << "Could not parse seed \"" << argv[i + 1] << "\" as integer.\n";
				err = 2;
				return false;
			}

			i++;
		} else if (!strcmp(argv[i], "-v")) {
			arg.verbose = true;
		} else {
			std::cout << "Unrecognised argument \"" << argv[i] << "\".\n";
			printHelp();
			err = 1;
			return false;
		}
	}

	return true;
}

void extractConfig(const json& configFile, const Arguments& arg, Configuration& config) {
	config.energyPerTransmit = configFile.at("energyPerTransmit");

	const json& nodes = configFile.at("nodes");
	config.nodes.resize(nodes.size());

	for (unsigned i = 0; i < nodes.size(); i++) {
		Node& node         = config.nodes[i];
		const json& jsNode = nodes.at(i);
		node.edgeNode      = jsNode.at("edge");
		node.energy = node.startingEnergy = jsNode.at("energy");
		node.observationRate              = jsNode.at("observationRate");
		node.dist                         = std::exponential_distribution<double>(node.observationRate);

		if (node.edgeNode) { config.edges.push_back(i); }
	}

	const json& links    = configFile.at("links");
	double defaultWeight = arg.alg == Algorithm::OSPF ? std::numeric_limits<double>::infinity() : 0;
	for (const json& link : links) {
		unsigned link1 = link.at(0);
		unsigned link2 = link.at(1);

		if (link1 >= config.nodes.size() || link2 >= config.nodes.size()) {
			throw std::out_of_range("Link references out of range node.");
		}

		config.nodes[link1].links.emplace_back(link2, config.nodes[link2].links.size(), defaultWeight);
		config.nodes[link2].links.emplace_back(link1, config.nodes[link1].links.size() - 1, defaultWeight);
	}
}

void printHelp() {
	Arguments arg;
	std::cout << "Usage: proj ospf <config> [options]                                          (1)\n"
	          << "Usage: proj ttl  <config> [options]                                          (2)\n"
	          << "   or: proj -h                                                               (3)\n\n"
	          << "(1) Run the project with the given configuration file using a standard OSPF\n"
	          << "    routing algorithm.\n"
	          << "(2) Run the project with the given configuration file using the proposed TTL\n"
	          << "    routing algorithm.\n"
	          << "(3) Print this help menu.\n\n"
	          << "OPTIONS\n"
	          << "  -s   <seed>  Set the seed used for generating random observations.\n"
	          << "  -v           Set the output to verbose. Prints initialization information\n"
	          << "               and final network state.\n";
}