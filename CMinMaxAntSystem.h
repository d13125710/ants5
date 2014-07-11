#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>

#include "CAnt.h"
#include "CAntSystem.h"
#include "CLocalSearch.h"




class CMinMaxAntSystem: public CAntSystem
{
private:
	double m_tau0;
	bool m_doLocalSearch , m_doOpt2 , m_doOpt3;
	//double m_minPheromone ,m_maxPheromone;
    int m_maxStagnationIterations , m_stagnationIterations;
	CLocalSearch *m_pLocalSearch;
	double lambda;
	double  m_branchingFactor;
	double  trail_0;
	double trail_max;
	double	trail_min;
	bool m_resetAnt;
	double m_restartBestAntTourLength;
	int m_bestDistanceLength;
	int m_bestIterationLength;
	int m_neg;

	std::vector<int> m_bestToDatePath;  //best to date for even
	std::vector<int> m_restartAntBestPath;
	std::vector<int> m_bestIterationAntPath;
	std::vector<int> m_iterationBestAntPath;

	//added
private:
	std::vector<std::vector<double> > m_nnList;
	void calculateNearestNeigbhor(int numberOfAnts);
	void initPheromoneTrails(double initialValue) const;
	void globalUpdatePheromone(const std::vector<int> &AntTour);
 	virtual void updateBestSoFarTour();
	void chooseClosestNext(std::vector<bool> &antsvisted , std::vector<int> &nntour);
	double nodeBranching(double l);

public:
	void localSearch();
	void initPheromones();
	void checkPheromoneLimits();
	void updatePheromones();
	CMinMaxAntSystem(int noNodes, int noAnts, int maxStagnationIterations , std::vector<std::vector<int> >  *Vdistance);
	virtual ~CMinMaxAntSystem(void);
};

