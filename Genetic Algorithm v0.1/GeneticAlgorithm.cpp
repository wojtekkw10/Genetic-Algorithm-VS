#include "GeneticAlgorithm.h"
#include <iostream>
double GeneticAlgorithm::CalculateAverageFitness()
{
	double sum = 0;
	for (int i = 0; i < Population.size(); i++)
	{
		sum += Population[i].fitness;
	}
	double average = sum / Population.size();
	avgFitness = average;
	return average;
}

std::vector <double> GeneticAlgorithm::GetTheBestChromosome()
{
	Sort();
	return Population[0].Chromosome;
}

void GeneticAlgorithm::Evolve(int Epochs)
{
	for (int i = 0; i < Epochs; i++)
	{
		epoch++;
		Sort();
		Mutate();

		CalculateFitness();
		CalculateAverageFitness();



		Sort();
		if (SelectionMethod == Roulette) RouletteSelection();
		if (SelectionMethod == Rank) RankSelection();
		Crossover();
	}



	if (displayStats)
	{
		CalculateFitness();
		Sort();
		std::cout << "Epoch: " << epoch << " ";
		std::cout << "The best fitness: " << Population[0].fitness << std::endl;
		std::cout << "Average Fitness: " << avgFitness << std::endl << std::endl;
	}
}

void GeneticAlgorithm::CalculateFitness()
{
	for (int i = 0; i < Population.size(); i++)
	{
		Population[i].fitness = FitnessFunc(Population[i].Chromosome);
	}
}

std::vector <Individual> GeneticAlgorithm::Sort()
{
	//sort the population by fitness
	std::vector <Individual> Sorted;
	std::vector <int> Checked; //Individuals ktore zostaly juz posortowane

	for (int k = 0; k < Population.size(); k++)
	{
		double theBestFitness = 99999999;
		int numberOfTheBest = 0;

		for (int i = 0; i < Population.size(); i++)
		{
			if (Population[i].fitness < theBestFitness)
			{
				bool isChecked = false;
				for (int j = 0; j < Checked.size(); j++)
				{
					if (Checked[j] == i) isChecked = true;
				}
				if (isChecked == false)
				{
					theBestFitness = Population[i].fitness;
					numberOfTheBest = i;
				}
			}

		}
		Checked.push_back(numberOfTheBest);
		Sorted.push_back(Population[numberOfTheBest]);
	}
	//Invert
	int j = 0;
	std::vector <Individual> QSorted;
	QSorted.resize(Sorted.size());
	for (int i = Sorted.size() - 1; i >= 0; i--)
	{
		QSorted[j] = Sorted[i];
		j++;
	}
	for (int i = 0; i < Sorted.size(); i++)
	{
		Sorted[i] = QSorted[i];
	}
	Population = Sorted;
	return Sorted;
	//Sorted
}

std::vector<Individual> GeneticAlgorithm::RouletteSelection()
{
	double FitnessSum = 0;
	for (int i = 0; i < Population.size(); i++)
	{
		FitnessSum += Population[i].fitness;
	}

	std::vector <Individual> Roulette;
	for (int i = 0; i < Population.size(); i++)
	{
		int r = rand() % (int)FitnessSum;
		int j = 0;
		double sum = 0;
		for (j = 0; j<Population.size()-1; j++)
		{
			sum += Population[j].fitness;
			if (sum > r) break;
		}
		Roulette.push_back(Population[j]);
	}
	Population = Roulette;
	return Roulette;
}

std::vector<Individual> GeneticAlgorithm::RankSelection()
{
	Sort();

	for (int i = 0, j=Population.size(); i < Population.size(); i++, j--)
	{
		Population[i].fitness = j;
	}
	RouletteSelection();

	return Population;
}

std::vector<Individual> GeneticAlgorithm::Crossover()
{
	std::vector <Individual> Crossover;
	//Dla kazdych dwoch
	for (int i = 0; i < Population.size() / 2; i++)
	{
		bool isCrossed = 0;
		if (rand() % 1000 < crossoverRate * 1000) isCrossed = true;
		else isCrossed = false;

		if (isCrossed)
		{
			for (int j = 0; j < 2; j++)
			{
				int crossing = rand() % Population[i].Chromosome.size(); //wylosuj punkt przeciecia chromosomu
				Individual temp;
				for (int k = 0; k < crossing; k++)
				{
					temp.Chromosome.push_back(Population[2 * i].Chromosome[k]);
				} //kopiuj geny do punktu crossing

				for (int k = crossing; k < Population[i].Chromosome.size(); k++)
				{
					temp.Chromosome.push_back(Population[2 * i + 1].Chromosome[k]);
				} //kopiuj geny od punktu crossing z kolejnego osobnika

				Crossover.push_back(temp);

			}
		}
		else
		{
			Crossover.push_back(Population[2 * i]);
			Crossover.push_back(Population[2 * i + 1]);
		}
		//Wygeneruj dwa nowe


	}
	Population = Crossover;
	return Crossover;
}

void GeneticAlgorithm::Mutate()
{
	for (int i = 0; i < Population.size(); i++)
	{
		if (rand() % 1000 < mutationRate * 1000) {
			//zamiana miejsc danych w obrebie chromosomu
			int random1 = rand() % Population[0].Chromosome.size();
			int random2 = rand() % Population[0].Chromosome.size();
			double temp = Population[i].Chromosome[random1];
			Population[i].Chromosome[random1] = Population[i].Chromosome[random2];
			Population[i].Chromosome[random2] = temp;

			//zmiana wartosci
			int random3 = rand() % Population[0].Chromosome.size();

			int random4 = rand() % Population.size();
			int random5 = rand() % Population[0].Chromosome.size();
			Population[i].Chromosome[random3] += Population[random4].Chromosome[random5]; //adding value from random chromosome in population
		}
	}
}

GeneticAlgorithm::GeneticAlgorithm(int sizeOfPopulation, int randomBottom, int randomUp, int sizeofChromosome)
{
	Population.resize(sizeOfPopulation);
	crossoverRate = 0.7;
	epoch = 0;
	avgFitness = 0;
	mutationRate = 0.02;
	displayStats = false;

	for (unsigned int i = 0; i < Population.size(); i++)
	{
		for (int j = 0; j < sizeofChromosome; j++)
		{
			Population[i].Chromosome.push_back((rand() % (randomUp-randomBottom)) + randomBottom);
		}
	}
}
GeneticAlgorithm::~GeneticAlgorithm()
{
}

Individual::Individual()
{
	fitness = 0;
}
Individual::~Individual()
{
}

PerceptronNeuron::PerceptronNeuron(int numberOfInputs, double threashold)
{
	Activation = 0;
	Output = 0;
	this->Threashold = threashold;
	this->numberOfInputs = numberOfInputs;
	for (int i = 0; i < numberOfInputs; i++)
	{
		Weights.push_back( ((double)rand() / (RAND_MAX)));//0 < r < 1
		//r = ((double) rand() / (RAND_MAX + 1)) -1 < r < 0?
	}

}
PerceptronNeuron::~PerceptronNeuron()
{
}

NeuronLayer::NeuronLayer( int numberOfNeurons)
{
	for (int i = 0; i < numberOfNeurons; i++)
	{
		PerceptronNeuron temp(numberOfNeurons);
		Layer.push_back(temp);
	}
}
NeuronLayer::~NeuronLayer()
{

}

void NeuralNetwork::SetInputs(std::vector<double> Inputs)
{
	this->Inputs.clear();
	this->Inputs.shrink_to_fit();
	for (int i = 0; i < Inputs.size(); i++)
	{
		this->Inputs.push_back(Inputs[i]);
	}
}

void NeuralNetwork::SetWeights(std::vector<double> Weights)
{
	int l = 0;
	for (int i = 0; i < NN.size(); i++)
	{
		for (int j = 0; j < NN[0].Layer.size(); j++)
		{
			for (int k = 0; k < NN[0].Layer[0].Weights.size(); k++)
			{
				NN[i].Layer[j].Weights[k] = Weights[l];
				l++;
			}
		}
	}
}

std::vector<double> NeuralNetwork::GetWeights()
{
	std::vector <double> Weights;
	for (int i = 0; i < NN.size(); i++)
	{
		for (int j = 0; j < NN[0].Layer.size(); j++)
		{
			for (int k = 0; k < NN[0].Layer[0].Weights.size(); k++)
			{
				Weights.push_back(NN[i].Layer[j].Weights[k]);
			}
		}
	}
	return Weights;
}

void NeuralNetwork::Evaluate()
{
	//i - number of weight in neuron
	//j = number of neuron in neuron layer
	//k - number of layer

	//Calculate activation values in first layer
	for (int j = 0; j < NN[0].Layer.size(); j++)
	{
		double Activation = 0; //sum of weights * inputs
		for (int i = 0; i < sizeofInput; i++)
		{
			Activation += Inputs[i] * NN[0].Layer[j].Weights[i]; // j - number of neuron in a layer
		}
		NN[0].Layer[j].Activation = Activation;
	}
	//Calculate Outputs (threashold, sigma function later)
	for (int j = 0; j < NN[0].Layer.size(); j++)
	{
		if (NN[0].Layer[j].Activation > NN[0].Layer[j].Threashold) NN[0].Layer[j].Output = 1;
		else NN[0].Layer[j].Output = 0;
	}
	//First Layer Done

	for (int k = 1; k < NN.size()-1; k++)
	{
		for (int j = 0; j < NN[k].Layer.size(); j++)
		{
			double Activation = 0; //sum of weights * inputs
			for (int i = 0; i < NN[k].Layer[j].numberOfInputs; i++)
			{
				Activation += NN[k-1].Layer[i].Output * NN[k].Layer[j].Weights[i];
			}
			NN[k].Layer[j].Activation = Activation;
		}
		for (int j = 0; j < NN[k].Layer.size(); j++)
		{
			if (NN[k].Layer[j].Activation > NN[k].Layer[j].Threashold) NN[k].Layer[j].Output = 1;
			else NN[k].Layer[j].Output = 0;
		}
	}

	//Evaluating the output layer

	int thelastLayer = NN.size()-1; //NN.size() is an output layer
	for (int j = 0; j < sizeofOutput; j++)
		{
			double Activation = 0; //sum of weights * inputs
			for (int i = 0; i < NN[thelastLayer].Layer[j].numberOfInputs; i++)
			{
				Activation += NN[thelastLayer - 1].Layer[i].Output * NN[thelastLayer].Layer[j].Weights[i];
			}
			NN[thelastLayer].Layer[j].Activation = Activation;
		}
	for (int j = 0; j < sizeofOutput; j++)
		{
			if (NN[thelastLayer].Layer[j].Activation > NN[thelastLayer].Layer[j].Threashold) NN[thelastLayer].Layer[j].Output = 1;
			else NN[thelastLayer].Layer[j].Output = 0;
		}


	for (int j = 0; j < sizeofOutput; j++)
	{
		Outputs.push_back(NN[thelastLayer].Layer[j].Output);
	}





	 
}

NeuralNetwork::NeuralNetwork(int numberOfLayers, int numberOfNeurons)
{
	for (int i = 0; i < numberOfLayers+1; i++) //+1 for output
	{
		NeuronLayer temp(numberOfNeurons);
		NN.push_back(temp);
	}
}
NeuralNetwork::~NeuralNetwork()
{

}