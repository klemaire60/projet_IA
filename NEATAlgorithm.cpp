#include "NEATAlgorithm.h"
#include <math.h>
#include "Simulation.h"
#include "SimulationJob.h"
#include <iostream>
#include "Environment.h"

#define WIN_SCORE 4
#define LOSE_SCORE -4
#define DRAW_SCORE 0.5


int NEATAlgorithm::population = 300;
double NEATAlgorithm::deltaDisjoint = 2.0;
double NEATAlgorithm::deltaWeights = 0.4;
double NEATAlgorithm::deltaThreshold = 1.0;

int NEATAlgorithm::staleSpecies = 15;

double NEATAlgorithm::mutateConnectionsChance = 0.25;
double NEATAlgorithm::perturbChance = 0.9;
double NEATAlgorithm::crossOverChance = 0.75;
double NEATAlgorithm::linkMutationChance = 2.0;
double NEATAlgorithm::nodeMutationChance = 0.5;
double NEATAlgorithm::biasMutationChance = 0.4;
double NEATAlgorithm::stepSize = 0.1;
double NEATAlgorithm::disableMutationChance = 0.4;
double NEATAlgorithm::enableMutationChance = 0.2;

double NEATAlgorithm::timeoutConstant = 20;

int NEATAlgorithm::maxNodes = 1000000;


NEATAlgorithm * NEATAlgorithm::instance = nullptr;

double NEATAlgorithm::sigmoid(double x)
{
	return 2.0 / (1 + exp(-4.9 * x)) - 1;
}

Genome * NEATAlgorithm::basicGenome()
{
	Genome * genome = new Genome();
	genome->maxNeuron = nbInputs;
	mutate(*genome);
	return genome;
}

void NEATAlgorithm::mutate(Genome & genome)
{
	for (std::map<std::string, double>::iterator it = genome.mutationRates.begin(); it != genome.mutationRates.end(); it++)
	{
		if (rand() % 2 == 1)
		{
			genome.mutationRates[(*it).first] = 0.95 * (*it).second;
		}
		else
		{
			genome.mutationRates[(*it).first] = 1.05263 * (*it).second;
		}
	}

	if (random() < genome.mutationRates["connections"])
	{
		pointMutate(genome);
	}

	double p = genome.mutationRates["link"];
	while (p > 0)
	{
		if (random() < p)
		{
			linkMutate(genome, false);
		}
		p -= 1;
	}

	p = genome.mutationRates["bias"];
	while (p > 0)
	{
		if (random() < p)
		{
			linkMutate(genome, true);
		}
		p -= 1;
	}

	p = genome.mutationRates["node"];
	while (p > 0)
	{
		if (random() < p)
		{
			nodeMutate(genome);
		}

		p -= 1;
	}

	p = genome.mutationRates["enable"];
	while (p > 0)
	{
		if (random() < p)
		{
			enableDisableMutate(genome, true);
		}

		p -= 1;
	}

	p = genome.mutationRates["disable"];
	while (p > 0)
	{
		if (random() < p)
		{
			enableDisableMutate(genome, false);
		}

		p -= 1;
	}
}

void NEATAlgorithm::pointMutate(Genome & genome)
{
	double step = genome.mutationRates["step"];

	for (int i = 0; i < genome.genes.size(); i++)
	{
		Gene * gene = genome.genes[i];

		if (random() < perturbChance)
		{
			gene->weight = gene->weight + random() * step * 2 - step;
		}
		else
			gene->weight = random() * 4 - 2;
	}
}

void NEATAlgorithm::linkMutate(Genome & genome, bool forceBias)
{
	int neuron1 = randomNeuron(genome.genes, false);
	int neuron2 = randomNeuron(genome.genes, true);

	Gene * newLink = new Gene();
	// Both are inputs :
	if (neuron1 < nbInputs && neuron2 < nbInputs)
	{
		return;
	}

	// Swap output and input :
	if (neuron2 < nbInputs)
	{
		int temp = neuron1;
		neuron1 = neuron2;
		neuron2 = temp;
	}

	newLink->into = neuron1;
	newLink->out = neuron2;

	if (forceBias)
	{
		newLink->into = nbInputs;
	}

	if (containsLink(genome.genes, *newLink))
	{
		return;
	}

	newLink->innovation = newInnovation();
	newLink->weight = random() * 4 - 2;
	genome.genes.push_back(newLink);
}

void NEATAlgorithm::nodeMutate(Genome & genome)
{
	if (genome.genes.size() == 0)
	{
		return;
	}

	genome.maxNeuron = genome.maxNeuron + 1;

	std::uniform_int_distribution<int> distribution(0, genome.genes.size() - 1);
	Gene * gene = genome.genes[distribution(device)];

	if (!gene->enabled)
		return;

	gene->enabled = false;

	Gene * gene1 = new Gene(*gene);
	gene1->out = genome.maxNeuron;
	gene1->weight = 1.0;
	gene1->innovation = newInnovation();
	gene1->enabled = true;
	
	Gene * gene2 = new Gene(*gene);
	gene2->into = genome.maxNeuron;
	gene2->innovation = newInnovation();
	gene2->enabled = true;

	genome.genes.push_back(gene1);
	genome.genes.push_back(gene2);
}

void NEATAlgorithm::enableDisableMutate(Genome & genome, bool enable)
{
	std::vector<Gene*> candidates;

	for (int i = 0; i < genome.genes.size(); i++)
	{
		if (genome.genes[i]->enabled != enable)
		{
			candidates.push_back(genome.genes[i]);
		}
	}

	if (candidates.size() == 0)
		return;

	std::uniform_int_distribution<int> distribution(0, candidates.size() - 1);
	Gene * gene = candidates[distribution(device)];
	gene->enabled = !gene->enabled;
}

double NEATAlgorithm::random()
{
	return (double)((double) rand() / (double)RAND_MAX);
}

bool NEATAlgorithm::containsLink(std::vector<Gene*>& genes, Gene & link)
{
	for (int i = 0; i < genes.size(); i++)
	{
		Gene * gene = genes[i];
		if (gene->into == link.into && gene->out == link.out) 
			return true;
	}

	return false;
}


// Take a random neuron in the genes :
int NEATAlgorithm::randomNeuron(std::vector<Gene*>& genes, bool nonInput)
{
	std::map<int, bool> neurons;
	if (!nonInput)
	{
		for (int i = 0; i < nbInputs; i++)
		{
			neurons[i] = true;
		}
	}

	for (int i = 0; i < nbOutputs; i++)
	{
		neurons[maxNodes + i] = true;
	}

	for (int i = 0; i < genes.size(); i++)
	{
		if (!nonInput || genes[i]->into > nbInputs)
		{
			neurons[genes[i]->into] = true;
		}

		if (!nonInput || genes[i]->out > nbInputs)
		{
			neurons[genes[i]->out] = true;
		}
	}

	std::uniform_int_distribution<int> distribution(0, neurons.size());
	int n = distribution(device);
	for (std::map<int, bool>::iterator it = neurons.begin(); it != neurons.end(); it++)
	{
		n--;
		if (n == 0)
		{
			return (*it).first;
		}
	}

	return 0;
}

int NEATAlgorithm::newInnovation()
{
	return ++innovation;
}

void NEATAlgorithm::generateNetwork(Genome & genome)
{
	Network * network = new Network();
	for(int i = 0; i < nbInputs; i++)
		network->neurons[i] = new NEATNeuron();

	for (int i = 0; i < nbOutputs; i++)
		network->neurons[maxNodes + i] = new NEATNeuron();

	std::sort(genome.genes.begin(), genome.genes.end(), [&](Gene * g1, Gene * g2) {
		return g1->out < g2->out;
	});

	for (int i = 0; i < genome.genes.size(); i++)
	{
		Gene * gene = genome.genes[i];
		if (network->neurons.find(gene->out) == network->neurons.end())
			network->neurons[gene->out] = new NEATNeuron();

		NEATNeuron * neuron = network->neurons[gene->out];
		neuron->incoming.push_back(gene);

		if (network->neurons.find(gene->into) == network->neurons.end())
			network->neurons[gene->into] = new NEATNeuron();
	}

	if (genome.network != nullptr) delete genome.network;
	genome.network = network;
}

void NEATAlgorithm::addToSpecies(Genome * child)
{
	bool foundSpecies = false;

	for (int i = 0; i < this->species.size(); i++)
	{
		Species * spec = species[i];
		if (!foundSpecies && sameSpecies(*spec->genomes[0], *child))
		{
			spec->genomes.push_back(child);
			foundSpecies = true;
			break;
		}
	}

	if (!foundSpecies)
	{
		Species * childSpecies = new Species();
		childSpecies->genomes.push_back(child);
		species.push_back(childSpecies);
	}
}

double NEATAlgorithm::disjoint(std::vector<Gene*>& genes1, std::vector<Gene*>& genes2)
{
	std::map<int, bool> innovation1;
	for (int i = 0; i < genes1.size(); i++)
	{
		Gene * gene = genes1[i];
		innovation1[gene->innovation] = true;
	}

	std::map<int, bool> innovation2;
	for (int i = 0; i < genes2.size(); i++)
	{
		Gene * gene = genes2[i];
		innovation2[gene->innovation] = true;
	}

	int disjointGenes = 0;

	for (int i = 0; i < genes1.size(); i++)
	{
		Gene * gene = genes1[i];
		if (innovation2.find(gene->innovation) == innovation2.end())
			disjointGenes++;
	}

	for (int i = 0; i < genes2.size(); i++)
	{
		Gene * gene = genes2[i];
		if (innovation1.find(gene->innovation) == innovation1.end())
			disjointGenes++;
	}

	double n = std::max<int>(genes1.size(), genes2.size());


	return (double)(((double) disjointGenes) / n);
}

double NEATAlgorithm::weights(std::vector<Gene*>& genes1, std::vector<Gene*>& genes2)
{
	std::map<int, Gene *> innovation2;
	for (int i = 0; i < genes2.size(); i++)
	{
		innovation2[genes2[i]->innovation] = genes2[i];
	}

	double sum = 0;
	double coincident = 0;

	for (int i = 0; i < genes1.size(); i++)
	{
		Gene * gene = genes1[i];

		if (innovation2.find(gene->innovation) != innovation2.end())
		{
			Gene * gene2 = innovation2[gene->innovation];
			sum += std::abs(gene->weight - gene2->weight);
			coincident += 1;
		}
	}

	if (coincident == 0)
		return 0;

	return sum / coincident;
}

bool NEATAlgorithm::sameSpecies(Genome & genome1, Genome & genome2)
{
	double dd = deltaDisjoint * disjoint(genome1.genes, genome2.genes);
	double dw = deltaWeights * weights(genome1.genes, genome2.genes);

	return dd + dw < deltaThreshold;
}

void NEATAlgorithm::cullSpecies(bool cutToOne)
{
	for (int i = 0; i < species.size(); i++)
	{
		Species * s = species[i];
		std::sort(s->genomes.begin(), s->genomes.end(), [&](Genome * g1, Genome * g2) {
			return g1->fitness > g2->fitness;
		});
		int remaining = std::ceil(s->genomes.size() / 2.0);

		if (cutToOne)
		{
			remaining = 1;
		}

		// TODO : Vérifier l'ordre du tri !
		while (s->genomes.size() > remaining)
		{
			Genome * last = s->genomes.back();
			s->genomes.pop_back();
			delete last;
		}
	}
}

void NEATAlgorithm::rankGlobally()
{
	std::vector<Genome*> global;

	for (int i = 0; i < species.size(); i++)
	{
		Species * spec = species[i];

		for (int j = 0; j < spec->genomes.size(); j++)
		{
			global.push_back(spec->genomes[j]);
		}
	}

	std::sort(global.begin(), global.end(), [&](Genome * g1, Genome * g2) {
		return g1->fitness < g2->fitness;
	});

	for (int i = 0; i < global.size(); i++)
	{
		global[i]->globalRank = i + 1;
	}
}

void NEATAlgorithm::removeStaleSpecies()
{
	std::vector<Species*> survived;

	for (int i = 0; i < species.size(); i++)
	{
		Species * s = species[i];

		std::sort(s->genomes.begin(), s->genomes.end(), [&](Genome * g1, Genome * g2) {
			return g1->fitness > g2->fitness;
		});

		if (s->genomes[0]->fitness > s->topFitness)
		{
			s->topFitness = s->genomes[0]->fitness;
			s->staleness = 0;
		}
		else
			s->staleness++;

		if (s->staleness < staleSpecies || s->topFitness >= maxFitness)
			survived.push_back(s);
		else
			delete s;
	}

	species = survived;
}

void NEATAlgorithm::removeWeakSpecies()
{
	std::vector<Species*> survived;

	double sum = totalAverageFitness();

	for (int i = 0; i < species.size(); i++)
	{
		Species * s = species[i];
		double breed = std::floor(s->averageFitness / sum * population);
		if (breed >= 1)
			survived.push_back(s);
		else
			delete s;
	}

	species = survived;
}

void NEATAlgorithm::calculateAverageFitness(Species & s)
{
	double total = 0;

	for (int i = 0; i < s.genomes.size(); i++)
	{
		total += s.genomes[i]->globalRank;
	}

	s.averageFitness = total / (double)s.genomes.size();
}

double NEATAlgorithm::totalAverageFitness()
{
	double total = 0;

	for (int i = 0; i < species.size(); i++)
	{
		total += species[i]->averageFitness;
	}

	return total;
}

void NEATAlgorithm::newGeneration()
{
	cullSpecies(false);
	rankGlobally();
	removeStaleSpecies();
	rankGlobally();
	for (int i = 0; i < species.size(); i++)
	{
		calculateAverageFitness(*species[i]);
	}
	removeWeakSpecies();
	double sum = totalAverageFitness();
	std::vector<Genome*> children;
	for (int i = 0; i < species.size(); i++)
	{
		Species * s = species[i];
		int breed = std::floor(s->averageFitness / sum * population) - 1;

		for (int i = 0; i < breed; i++)
		{
			children.push_back(breedChild(*s));
		}
	}
	cullSpecies(true);

	std::uniform_int_distribution<int> distribution(0, species.size() - 1);
	while (children.size() + species.size() < population)
	{
		Species * s = species[distribution(device)];
		children.push_back(breedChild(*s));
	}
	for (int i = 0; i < children.size(); i++)
	{
		addToSpecies(children[i]);
	}
	std::vector<Genome*> genomes = getAllGenomes();
	for (int i = 0; i < genomes.size(); i++)
	{
		generateNetwork(*genomes[i]);
	}
	generation++;
}

Genome * NEATAlgorithm::breedChild(Species & s)
{
	Genome * child;

	std::uniform_int_distribution<int> distribution(0, s.genomes.size() - 1);

	if (random() < crossOverChance)
	{
		//std::cout << "Crossover" << std::endl;
		Genome * g1 = s.genomes[distribution(device)];
		Genome * g2 = s.genomes[distribution(device)];
		child = crossover(g1, g2);
	}
	else
	{
		//std::cout << "Not crossover" << std::endl;
		child = new Genome(*s.genomes[distribution(device)]);
	}

	generateNetwork(*child);
	mutate(*child);
	generateNetwork(*child);

	return child;
}

Genome * NEATAlgorithm::crossover(Genome * g1, Genome * g2)
{
	Genome * result = new Genome();

	if (g2->fitness > g1->fitness)
	{
		Genome * tmp = g1;
		g1 = g2;
		g2 = tmp;
	}

	std::map<int, Gene*> innovations2;
	for (int i = 0; i < g2->genes.size(); i++)
	{
		Gene * gene = g2->genes[i];
		innovations2[gene->innovation] = gene;
	}

	for (int i = 0; i < g1->genes.size(); i++)
	{
		Gene * gene1 = g1->genes[i];

		if (innovations2.find(gene1->innovation) != innovations2.end())
		{
			Gene * gene2 = innovations2[gene1->innovation];
			result->genes.push_back(new Gene(*gene2));
		}
		else
		{
			result->genes.push_back(new Gene(*gene1));
		}
	}

	result->maxNeuron = std::max<int>(g1->maxNeuron, g2->maxNeuron);

	for (std::map<std::string, double>::iterator it = g1->mutationRates.begin(); it != g1->mutationRates.end(); it++)
	{
		result->mutationRates[(*it).first] = (*it).second;
	}

	return result;
}

void NEATAlgorithm::initializePool()
{
	for (int i = 0; i < population; i++)
	{
		Genome * child = basicGenome();
		generateNetwork(*child);
		mutate(*child);
		generateNetwork(*child);
		addToSpecies(child);
	}
}

std::vector<double> NEATAlgorithm::evaluateNetwork(Network & network, std::vector<double> inputs)
{
	std::vector<double> result;

	network.mutex.lock();
	for (int i = 0; i < inputs.size(); i++)
	{
		network.neurons[i]->value = inputs[i];
	}

	for (std::map<int, NEATNeuron*>::iterator it = network.neurons.begin(); it != network.neurons.end(); it++)
	{
		double sum = 0;
		NEATNeuron * neuron = (*it).second;
		for (int j = 0; j < neuron->incoming.size(); j++)
		{
			Gene * incoming = neuron->incoming[j];
			NEATNeuron * other = network.neurons[incoming->into];
			sum += incoming->weight * other->value;
		}

		if (neuron->incoming.size() > 0)
			neuron->value = sigmoid(sum);
	}

	for (int i = 0; i < NEATAlgorithm::getInstance()->nbOutputs; i++)
	{
		result.push_back(network.neurons[maxNodes + i]->value);
	}
	network.mutex.unlock();

	return result;
}

NEATAlgorithm * NEATAlgorithm::getInstance()
{
	return instance;
}

Genome * NEATAlgorithm::getNeuralNetwork(int nnid)
{
	return bestGenome;
}


std::vector<Genome*> NEATAlgorithm::getAllGenomes()
{
	std::vector<Genome*> generation;
	for (int i = 0; i < species.size(); i++)
	{
		for (int j = 0; j < species[i]->genomes.size(); j++)
		{
			generation.push_back(species[i]->genomes[j]);
		}
	}
	return generation;
}

/*
void NEATAlgorithm::check()
{
	std::vector<Genome*> gens = getAllGenomes();
	bool hasSameAddress = false;
	for (int i = 0; i < gens.size(); i++)
	{
		for (int j = i + 1; j < gens.size(); j++)
		{
			if (gens[i] == gens[j] || gens[i]->network == gens[j]->network)
			{
				hasSameAddress = true;
				break;
			}
		}
	}

	for (int i = 0; i < gens.size(); i++)
	{
		std::cout << gens[i]->network->neurons.size() << std::endl;
		if (gens[i]->network->neurons.size() > 100)
		{
			std::cout << gens[i]->network->neurons[0]->value << std::endl;
		}
	}
}
*/

NEATAlgorithm::NEATAlgorithm(int nbInputs, int nbOutputs)
{
	/*
	population = 300;
	deltaDisjoint = 2.0;
	deltaWeights = 0.4;
	deltaThreshold = 1.0;

	staleSpecies = 15;

	mutateConnectionsChance = 0.25;
	perturbChance = 0.9;
	crossOverChance = 0.75;
	linkMutationChance = 2.0;
	nodeMutationChance = 0.5;
	biasMutationChance = 0.4;
	stepSize = 0.1;
	disableMutationChance = 0.4;
	enableMutationChance = 0.2;

	timeoutConstant = 20;

	maxNodes = 1000000;
	*/

	this->nbInputs = nbInputs;
	this->nbOutputs = nbOutputs;

	generation = 0;
	innovation = nbOutputs;	// nb outputs
	currentSpecies = 1;
	currentGenome = 1;
	currentFrame = 0;
	maxFitness = -1000000;
	instance = this;

	bestGenome = nullptr;

	initializePool();
}

void NEATAlgorithm::evaluatePopulation()
{
	std::vector<Genome*> population = getAllGenomes();
	for (int i = 0; i < population.size(); i++)
	{
		population[i]->fitness = 0;
	}
	
	std::cout << "Simulate game ..." << std::endl;
	// Create environment ...
	Environment env;
	for (int i = 0; i < population.size(); i++)
	{
		Robot * actor = new Robot();
		actor->brain = population[i];
		env.addActor(actor);
	}
	
	while (env.haveActorAlive())
	{
		env.iterate();
		std::vector<Robot *> & actors = env.getActors();
		double currentX = actors[0]->x;
		if (currentX > 10000)
		{
			for (int i = 0; i < actors.size(); i++)
			{
				if (!actors[i]->isDead)
				{
					actors[i]->fitness = currentX;
				}
			}
			break;
		}
	}

	std::vector<Robot*> actors = env.getActors();
	for (int i = 0; i < actors.size(); i++)
	{
		actors[i]->brain->fitness = actors[i]->fitness;
		delete actors[i];
	}

	printf("Simulation completed, check evaluations ...\n");


	std::sort(population.begin(), population.end(), [&](Genome *p1, Genome *p2) {
		return p1->fitness > p2->fitness;
	});

	bestGenome = population[0];

	for (int i = 0; i < population.size(); i++)
	{
		if (population[i]->fitness > 100)
		{
			std::cout << "NN" << i << " fitness : " << population[i]->fitness << std::endl;
		}
	}

	std::cout << "Best fitness : " << bestGenome->fitness << std::endl;
	std::cout << "Nb neuron : " << bestGenome->network->neurons.size() - 3 - 1 << std::endl;
}

Species::Species()
{
	topFitness = -1000000;
	staleness = 0;
	averageFitness = -1000000;
}

Species::~Species()
{
	for (int i = 0; i < genomes.size(); i++)
	{
		delete genomes[i];
	}
}

Genome::Genome()
{
	fitness = -1000000;
	adjustedFitness = -1000000;
	maxNeuron = 0;
	globalRank = 0;
	network = nullptr;
	mutationRates["connections"] = NEATAlgorithm::mutateConnectionsChance;
	mutationRates["link"] = NEATAlgorithm::linkMutationChance;
	mutationRates["bias"] = NEATAlgorithm::biasMutationChance;
	mutationRates["node"] = NEATAlgorithm::nodeMutationChance;
	mutationRates["enable"] = NEATAlgorithm::enableMutationChance;
	mutationRates["disable"] = NEATAlgorithm::disableMutationChance;
	mutationRates["step"] = NEATAlgorithm::stepSize;
}

Genome::Genome(const Genome & toCopy)
{
	for (int i = 0; i < toCopy.genes.size(); i++)
	{
		genes.push_back(new Gene(*toCopy.genes[i]));
	}

	maxNeuron = toCopy.maxNeuron;
	this->fitness = toCopy.fitness;
	this->adjustedFitness = toCopy.adjustedFitness;
	this->globalRank = toCopy.globalRank;

	this->network = nullptr;

	for (std::map<std::string, double>::const_iterator it = toCopy.mutationRates.cbegin(); it != toCopy.mutationRates.cend(); it++)
	{
		mutationRates[(*it).first] = (*it).second;
	}
}

Genome::~Genome()
{
	for (int i = 0; i < genes.size(); i++)
	{
		delete genes[i];
	}

	if(network != nullptr)
		delete network;
}

std::vector<double> Genome::evaluate(std::vector<double> inputs)
{
	return NEATAlgorithm::evaluateNetwork(*network, inputs);
}

void Genome::addFitness(double value)
{
	this->fitness += value;
}

Gene::Gene()
{
	into = 0;
	out = 0;
	weight = 0;
	enabled = true;
	innovation = 0;
}

Gene::Gene(const Gene & toCopy)
{
	into = toCopy.into;
	out = toCopy.out;
	weight = toCopy.weight;
	enabled = toCopy.enabled;
	innovation = toCopy.innovation;
}

NEATNeuron::NEATNeuron()
{
	value = 0;
}

NEATNeuron::NEATNeuron(const NEATNeuron & toCopy)
{
	for (int i = 0; i < toCopy.incoming.size(); i++)
		incoming.push_back(toCopy.incoming[i]);
}

NEATNeuron::~NEATNeuron()
{
}

Network::Network()
{
}

/*
Network::Network(const Network & toCopy)
{
	for (std::map<int, NEATNeuron*>::const_iterator it = toCopy.neurons.cbegin(); it != toCopy.neurons.cend(); it++)
	{
		neurons[(*it).first] = new NEATNeuron(*(*it).second);
	}
}
*/

Network::~Network()
{
	for (std::map<int, NEATNeuron*>::iterator it = neurons.begin(); it != neurons.end(); it++)
	{
		delete (*it).second;
	}
}

/*
void Network::operator=(const Network & toCopy)
{
	for (std::map<int, NEATNeuron*>::const_iterator it = toCopy.neurons.cbegin(); it != toCopy.neurons.cend(); it++)
	{
		neurons[(*it).first] = new NEATNeuron(*(*it).second);
	}
}
*/