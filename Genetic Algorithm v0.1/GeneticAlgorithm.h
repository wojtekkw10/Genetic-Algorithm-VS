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
	GeneticAlgorithm(int sizeOfPopulation, int randomBottom, int randomUp);
	~GeneticAlgorithm();
};

class Neuron {
	double weight;

	Neuron();
	~Neuron();
};