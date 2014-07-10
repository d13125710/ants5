#pragma once

#include "Chromo.h"
#include "LocalSearch.h"

typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class GeneticSystem
{
	int populationSize , noNodes;
    double mutationRatio;
    bool doOpt2;
    bool doOpt3;
	std::vector<Chromo> population;
    MatrixArrayTypeInt *dist;
    std::vector<double>  fitness;
    std::vector<int>  bestSoFarTour;
    int iteration;
    LocalSearch *opt;
	int firsttime;

	void updateBestSoFarTour();
	void computeFitness();
	void mutatePopulation(std::vector<Chromo> &pop);
	void sortPopulationAscending(std::vector<Chromo> &pop);
	void sortPopulationDescending(std::vector<Chromo> &pop);
	int tournamentSelection();
	void stepGeneration2();
	std::vector<int> &getBestTour();
	int getBestTourLength();
	int getBestSoFarTourLength();


public:

	GeneticSystem(int populationSize, int noNodes, double mutationRatio, bool doOpt2, bool doOpt3 , MatrixArrayTypeInt *d);
  	void setIteration(int iter);
	int getIteration();
	void initPopulation();
	void stepGeneration();
 	std::vector<int>  &getBestSoFarTour();
	void localSearch();
	int computeTourLength(std::vector<int> &tour);
	GeneticSystem(void);
	~GeneticSystem(void);


private:
	 int getRandomNumber(int low, int high, int except)
	 {
		//Generate random value between low and high (inclusive), except 'except'.
		int temp;
		bool done = false;

		while(!done)
		{
			temp = (rand() % ((high + 1) - low)) + low;
			if(temp != except)
			{
				done = true;
			}
		}
		return temp;
	 }


};

