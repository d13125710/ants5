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
	virtual void updateBestSoFarTour();
	virtual void updatePheromones();
	virtual void calculateHeuristicMatrix();
	virtual void decisionRule(int k, int step);
	virtual void evaporateAllPhero();
	
	void globalPheromoneDeposit();
	void localPheromoneUpdate(int ant, int step);
	virtual void initPheromones();
	CAntColonySystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix);
	~CAntColonySystem(void);
};

