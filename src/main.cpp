#include "main.h"

int main(int argc, char** argv) {
	int err;
	Arguments arg;

	if (!verifyArguments(argc, argv, arg, err)) { return err; }

	Configuration config;
	try {
		extractConfig(arg.config, arg, config);
	} catch (json::exception e) {
		std::cout << "Given configuration file was invalid:\n" << e.what() << std::endl;
		return 2;
	}

	Output out;
	switch (arg.alg) {
		case Algorithm::OSPF:
			simOSPF(arg, config, out);
			break;
		case Algorithm::TTL:
			simTTL(arg, config, out);
			break;
	}

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

	for (Node& node : config.nodes) {
		node.shortestLink = std::distance(
		    node.links.begin(), std::max_element(node.links.begin(), node.links.end(), std::greater<Connection>()));
	}

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

	for (unsigned i = 0; i < config.nodes.size(); i++) {
		Node& node = config.nodes[i];
		observationEvents.emplace_back(node.dist(gen), i);
	}

	std::make_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

	if (arg.verbose) {
		std::cout << "Initial Observations: (heap)\n"
		          << "Time       Node\n";
		for (const std::pair<double, unsigned>& observation : observationEvents) {
			printf("%8f   %4d\n", observation.first, observation.second);
		}
	}

	out.time                = 0;
	out.escapedObservations = 0;
	while (true) {
		std::pair<double, unsigned> observation = observationEvents.front();
		Node& observingNode                     = config.nodes[observation.second];
		std::pop_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

		double deltaT = observation.first;
		out.time += deltaT;
		std::for_each(observationEvents.begin(), observationEvents.end(),
		              [deltaT](std::pair<double, unsigned>& observation) { observation.first -= deltaT; });

		observationEvents.back() = {observingNode.dist(gen), observation.second};
		std::push_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());
		observingNode.observationsMade++;

		unsigned currentNodeIndex = observation.second;
		while (!config.nodes[currentNodeIndex].edgeNode &&
		       config.nodes[currentNodeIndex].energy >= config.energyPerTransmit) {
			Node& currentNode = config.nodes[currentNodeIndex];

			currentNode.energy -= config.energyPerTransmit;
			currentNodeIndex = currentNode.links[currentNode.shortestLink].nodeIndex;
		}

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

			double ttl1 = energy / config.energyPerTransmit / (observationRate + node.observationRate);
			double ttl2 = node.energy / config.energyPerTransmit / node.observationRate;
			double ttl  = std::min(ttl1, ttl2);

			Connection& con = node.links[connectionIndex];
			if (con.weight >= ttl) { continue; }

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

	for (Node& node : config.nodes) {
		node.shortestLink = std::distance(node.links.begin(), std::max_element(node.links.begin(), node.links.end()));
	}

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

	for (unsigned i = 0; i < config.nodes.size(); i++) {
		Node& node = config.nodes[i];
		observationEvents.emplace_back(node.dist(gen), i);
	}

	std::make_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

	if (arg.verbose) {
		std::cout << "Initial Observations: (heap)\n"
		          << "Time       Node\n";
		for (const std::pair<double, unsigned>& observation : observationEvents) {
			printf("%8f   %4d\n", observation.first, observation.second);
		}
	}

	out.time                = 0;
	out.escapedObservations = 0;
	while (true) {
		std::pair<double, unsigned> observation = observationEvents.front();
		Node& observingNode                     = config.nodes[observation.second];
		std::pop_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());

		double deltaT = observation.first;
		out.time += deltaT;
		std::for_each(observationEvents.begin(), observationEvents.end(),
		              [deltaT](std::pair<double, unsigned>& observation) { observation.first -= deltaT; });

		observationEvents.back() = {observingNode.dist(gen), observation.second};
		std::push_heap(observationEvents.begin(), observationEvents.end(), std::greater<std::pair<double, unsigned>>());
		observingNode.observationsMade++;

		unsigned currentNodeIndex = observation.second;
		std::stack<unsigned> routedPath;
		while (!config.nodes[currentNodeIndex].edgeNode &&
		       config.nodes[currentNodeIndex].energy >= config.energyPerTransmit) {
			Node& currentNode = config.nodes[currentNodeIndex];

			currentNode.energy -= config.energyPerTransmit;
			routedPath.push(currentNodeIndex);
			currentNodeIndex = currentNode.links[currentNode.shortestLink].nodeIndex;

			double deltaTSinceLastRoute =
			    config.nodes[currentNodeIndex].routedObservations.size() > 0 ?
                    out.time - config.nodes[currentNodeIndex].routedObservations.back().first :
                    out.time;

			config.nodes[currentNodeIndex].routedObservations.emplace_back(out.time, deltaTSinceLastRoute);
		}

		Node& lastNode = config.nodes[currentNodeIndex];
		if (lastNode.edgeNode && lastNode.energy >= config.energyPerTransmit) {
			lastNode.energy -= config.energyPerTransmit;
			out.escapedObservations++;

			double lastEnergy;
			double lastRate;

			while (!routedPath.empty()) {
				Node& currentNode = config.nodes[currentNodeIndex];
				currentNodeIndex  = routedPath.top();
				routedPath.pop();
				Node& prevNode = config.nodes[currentNodeIndex];

				double perceivedRouteReceiveRate =
				    std::distance(
				        currentNode.routedObservations.crbegin(),
				        std::find_if(
				            currentNode.routedObservations.crbegin(), currentNode.routedObservations.crend(),
				            [&out](const std::pair<double, int>& obs) { return out.time - obs.first < TTL_WINDOW; })) /
				    std::min(out.time, TTL_WINDOW);

				double ttl1 = currentNode.energy / config.energyPerTransmit /
				              (perceivedRouteReceiveRate + currentNode.observationRate);
				double ttl2 = currentNode.edgeNode ?
                                  std::numeric_limits<double>::infinity() :
                                  lastEnergy / config.energyPerTransmit /
				                      (perceivedRouteReceiveRate + currentNode.observationRate + lastRate);

				double deltaTSinceLastRoute = currentNode.routedObservations.back().second;

				std::for_each(
				    prevNode.links.begin(), prevNode.links.end(), [deltaTSinceLastRoute, &out](Connection& con) {
					    con.weight -= std::min(
					        con.weight, pow(TTL_BIAS, out.time - con.lastTimeWeightUpdated) * deltaTSinceLastRoute);
				    });

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