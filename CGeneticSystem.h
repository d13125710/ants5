#pragma once

#include "CChromo.h"
#include "CLocalSearch.h"

typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class CGeneticSystem
{
	int m_populationSize , m_noCitys;
    double m_mutationRatio;
    bool m_bOpt2Search;;
    bool m_bOpt3Search;
	std::vector<CChromo> m_ChromoPopulation;
    MatrixArrayTypeInt *m_distMatrix;
    std::vector<double>  m_fitness;
    std::vector<int>  m_bestSoFarPath;
    int m_currentIteration;
    CLocalSearch *m_pLocalSearch;
	int firsttime;

	void updateBestSoFarPath();
	void computeFitness();
	void mutatePopulation(std::vector<CChromo> &pop);
	void SortPopulation(std::vector<CChromo> &pop , bool Direction=true);

	int tournamentSelection();
	void stepGeneration2();
	std::vector<int> &getBestPath();
	int getBestPathLength();
	int getBestSoFarPathLength();


public:

	CGeneticSystem(int populationSize, int noNodes, double mutationRatio, bool doOpt2, bool doOpt3 , MatrixArrayTypeInt *d);
  	void setIteration(int iter);
	int getIteration() const;
	void initPopulation();
	void stepGeneration();
 	const std::vector<int> &getBestSoFarPath();
	void localSearch();
	int computePathLength(const std::vector<int> &tour);
	CGeneticSystem(void);
	~CGeneticSystem(void);


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

