#include <iostream>
#include "GeneticAlgorithm.h"
#include <time.h>


double CalculateFitness(std::vector <double> Chromosome)
{
	double sum = 0;
	for (int j = 0; j < Chromosome.size(); j++)
	{
		sum += Chromosome[j];
	}
	if (250 == sum) return 2000;
	else return -abs(250 - sum) + 2000;
}

int main()
{
	srand((unsigned int)time(NULL));

	for (int k = 0; k < 1; k++)
	{
		GeneticAlgorithm alg(100, 0, 500);
		alg.crossoverRate = 0.8;
		alg.mutationRate = 0.05;
		alg.FitnessFunc = CalculateFitness;
		alg.displayStats = true;
		alg.SelectionMethod = Rank;

		bool stop = false;
		for (int i = 0; i < 5; i++)
		{
			alg.Evolve(5);
			if (alg.Population[0].fitness == 2000) break;
		}


		//wypisz najlepszego
		std::cout << k << ". ";
		double sum = 0;
		for (int i = 0; i < 5; i++)
		{
			sum += alg.Population[0].Chromosome[i];
			std::cout << alg.Population[0].Chromosome[i] << " ";
		}
		std::cout << sum << " ";
		std::cout << std::endl << std::endl;
	}

	NeuralNetwork NN(6,20);
	std::cout << "fefefe";
	NN.Evaluate();
	/*
	alg.evolve();//crossover
	
	neuralnetwork neuralnetwork(alg.pop[n].chromosome//weights);
	alg.pop[n].chromosome = neuralnetwork.backpropagate(); //dodatkowo polepszamy

	alg.Population[0].Chromosome = neuralnetwork.fitness(alg.pop[n].chromosome//weights); //fitness in NN is 1/error
	alg.calculate;//crossover
	//now we have a population of better weights

	//NN przechowuje jeden chromosome(wagi) potrafi zrobic backpropagation i wyliczyc error/fitness
	//GA przechowuje wszystkie chromosomy populacji, potrafi je crossoverowac, mutowac itd
	*/
	system("pause");
}

//GA ktora generuje liczby ktorych suma odpowiada wartosci podanej przez usera

//Implement crossover rate

//Przesylanie funkcji obliczajacej fitness do algorutmu genetycznego zeby mogl go wykozystac np przed koncowym sortowaniem

//Multiple population in single GA