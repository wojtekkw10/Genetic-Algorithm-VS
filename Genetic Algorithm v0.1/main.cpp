#include <iostream>
#include "GeneticAlgorithm.h"
#include <time.h>

NeuralNetwork NN(1,4);

double CalculateFitnessNN(std::vector <double> Chromosome)
{
	/* Goal - for numbers less than 0 - 0; otherwise 1*/
	std::vector <double> Input;
	Input.push_back(0.5);

	NN.sizeofOutput = 1;
	NN.sizeofInput = 1;
	NN.SetWeights(Chromosome);
	NN.SetInputs(Input);
	NN.Evaluate();
	std::vector <double> Output = NN.Outputs;

	double error = 0;

	error += pow(Output[0] - 1, 2);

	if (0 == error) return 10;
	else return -abs(error) + 10;

}


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
		GeneticAlgorithm alg(100, 0, 2, 100);
		alg.crossoverRate = 0.8;
		alg.mutationRate = 0.05;
		alg.FitnessFunc = CalculateFitnessNN;
		alg.displayStats = true;
		alg.SelectionMethod = Rank;

		bool stop = false;
		for (int i = 0; i < 5; i++)
		{
			alg.Evolve(1);
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
	std::cout << NN.Outputs[0];


	std::cout << std::endl << std::endl << std::endl << std::endl;

	for (int i = 0; i < 100; i++)
	{
		std::vector <double> Input;
		Input.push_back(0.51);

		NeuralNetwork NNN(5, 20);
		NNN.sizeofOutput = 1;
		NNN.sizeofInput = 1;
		NNN.SetInputs(Input);
		NNN.Evaluate();
		std::cout << NNN.Outputs[0] << std::endl;
	}
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