#pragma once
#include <vector>
#include <map>
#include <string>
#include <random>
#include <mutex>

struct Gene
{
	int into;
	int out;
	double weight;
	bool enabled;
	int innovation;

	Gene();
	Gene(const Gene & toCopy);
};

struct NEATNeuron
{
	std::vector<Gene*> incoming;
	double value;

	NEATNeuron();
	NEATNeuron(const NEATNeuron & toCopy);
	~NEATNeuron();
};

struct Network
{
	std::map<int, NEATNeuron*> neurons;
	std::mutex mutex;

	Network();
	//Network(const Network & toCopy);
	~Network();
	//void operator=(const Network & toCopy);
};

struct Genome
{
	std::vector<Gene*> genes;
	double fitness;
	double adjustedFitness;
	Network * network;
	int maxNeuron;
	int globalRank;
	std::map<std::string, double> mutationRates;

	Genome();
	Genome(const Genome & toCopy);
	~Genome();

	virtual std::vector<double> evaluate(std::vector<double> inputs);
	virtual void addFitness(double value);

};

struct Species
{
	double topFitness;
	double staleness;
	std::vector<Genome*> genomes;
	double averageFitness;

	Species();
	~Species();
};

class NEATAlgorithm
{
	static NEATAlgorithm * instance;

	int nbInputs;
	int nbOutputs;

	std::vector<Species*> species;
	int generation;
	int innovation;
	int currentSpecies;
	int currentGenome;
	int currentFrame;
	double maxFitness;

	std::random_device device;

	static double sigmoid(double x);
	Genome * basicGenome();
	void mutate(Genome & genome);
	void pointMutate(Genome & genome);
	void linkMutate(Genome & genome, bool forceBias);
	void nodeMutate(Genome & genome);
	void enableDisableMutate(Genome & genome, bool enable);

	double random();
	bool containsLink(std::vector<Gene*> & genes, Gene & link);
	int randomNeuron(std::vector<Gene*> & genes, bool nonInput);

	int newInnovation();

	void generateNetwork(Genome & genome);
	void addToSpecies(Genome * child);

	double disjoint(std::vector<Gene*> & genes1, std::vector<Gene*> & genes2);
	double weights(std::vector<Gene*> & genes1, std::vector<Gene*> & genes2);
	bool sameSpecies(Genome & genome1, Genome & genome2);

	void cullSpecies(bool cutToOne);
	void rankGlobally();
	void removeStaleSpecies();
	void removeWeakSpecies();
	void calculateAverageFitness(Species & s);
	double totalAverageFitness();
	Genome * breedChild(Species & s);
	Genome * crossover(Genome * g1, Genome * g2);
	void initializePool();
	
	Genome * bestGenome;

	//void check();

public:
	std::vector<Genome*> getAllGenomes();
	NEATAlgorithm(int nbInputs, int nbOutputs);
	void evaluatePopulation();
	static std::vector<double> evaluateNetwork(Network & network, std::vector<double> inputs);
	static NEATAlgorithm * getInstance();
	Genome * getNeuralNetwork(int nnid);

	void newGeneration();
	int getGenerationNumber() { return generation; }


	static int population;
	static double deltaDisjoint;
	static double deltaWeights;
	static double deltaThreshold;

	static int staleSpecies;

	static double mutateConnectionsChance;
	static double perturbChance;
	static double crossOverChance;
	static double linkMutationChance;
	static double nodeMutationChance;
	static double biasMutationChance;
	static double stepSize;
	static double disableMutationChance;
	static double enableMutationChance;

	static double timeoutConstant;

	static int maxNodes;
};

