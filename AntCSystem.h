#pragma once
#include "antsystem.h"
#include "LocalSearch.h"

class AntCSystem :
	public AntSystem
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


	virtual void Evaporation();
	void globalPheromoneDeposit();
	void localPheromoneUpdate(int ant, int step);
	virtual void initPheromones();
	double fRand(double fMin, double fMax)
	{
		double f = (double)rand() / RAND_MAX;
		return fMin + f * (fMax - fMin);
	}


	AntCSystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix);
	~AntCSystem(void);
};

