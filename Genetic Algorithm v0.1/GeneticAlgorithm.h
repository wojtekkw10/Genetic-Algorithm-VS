#include <vector>
#include <functional>
#pragma once
enum SelectionMethod
{
	Roulette,
	Rank
};


class Individual
{
private:

public:
	std::vector <double> Chromosome;
	double fitness;

	Individual();
	~Individual();
};

class GeneticAlgorithm
{
private:
	
public:
	std::vector <Individual> Population;
	int numberOfPopulations;
	double crossoverRate;
	double mutationRate;
	int epoch;
	double avgFitness;
	bool displayStats;
	std::function<double(std::vector <double> Chromosome)> FitnessFunc; //Wskaznik do funkcji, ktora przyjmuje chromosom i zwraca fitness

	void Evolve(int Epochs = 1);
	void CalculateFitness();
	void Mutate();
	double CalculateAverageFitness();
	std::vector <double> GetTheBestChromosome();

	std::vector <Individual> Sort();
	std::vector <Individual> RouletteSelection();
	std::vector <Individual> RankSelection();
	SelectionMethod SelectionMethod;

	std::vector <Individual> Crossover();
	GeneticAlgorithm(int sizeOfPopulation, int randomBottom, int randomUp, int sizeofChromosome);
	~GeneticAlgorithm();
};

class PerceptronNeuron {
public:
	std::vector <double> Weights;
	int numberOfInputs;
	double Threashold;
	double Activation;
	double Output;

	PerceptronNeuron(int numberOfInputs, double threashold = 0.5);
	~PerceptronNeuron();
};

class NeuronLayer
{
public:
	std::vector <PerceptronNeuron> Layer;

	NeuronLayer(int numberOfNeurons);
	~NeuronLayer();
};

class NeuralNetwork
{
public:
	std::vector <NeuronLayer> NN;
	std::vector <double> Inputs; //Place for inputs
	std::vector <double> Outputs;

	int sizeofOutput;
	int sizeofInput;

	void SetInputs(std::vector <double> Inputs); //Get inputs for evaluation

	void SetWeights(std::vector <double> Weights); //from Genetic algorithm as an array of weight
	std::vector <double> GetWeights(); //for GA

	void Evaluate(); //Calculate the outputs from Inputs and weights

	NeuralNetwork(int numberOfLayers, int numberOfNeurons);
	~NeuralNetwork();
};