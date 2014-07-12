#pragma once
#include "CAntSystem.h"
#include "CLocalSearch.h"
#include "Utils.h"

class CAntColonySystem :
	public CAntSystem
{
	  double xi;
      double tau0;
      double q0;
	  double m_r0;
	

public:

	virtual void constructSolutions();
	virtual void updateBestSoFarPath();
	virtual void updatePheromones();
	virtual void calculateHeuristicMatrix();
	virtual void decisionRule(int k, int step);
	virtual void evaporateAllPhero();

	double HEURISTIC(int m, int n) {
		return (1.0 / ((double) (*m_distanceMatrix)[m][n] + 0.1));
	}
	
	void globalPheromoneDeposit();
	void localPheromoneUpdate(int ant, int step);
	virtual void initPheromones();
	CAntColonySystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix);
	~CAntColonySystem(void);


	void calculateNearestNeigbhor(unsigned int numberOfAnts);
	std::vector<std::vector<int> > m_nnList;
	void chooseClosestNext(std::vector<bool> &antsvisted , std::vector<int> &nntour);
};

