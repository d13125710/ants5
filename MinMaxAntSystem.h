#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>

#include "Ant.h"
#include "AntSystem.h"
#include "LocalSearch.h"




class MinMaxAntSystem: public AntSystem
{
	double m_tau0;
	bool m_doLocalSearch , m_doOpt2 , m_doOpt3;
	//double m_minPheromone ,m_maxPheromone;
    int m_maxStagnationIterations , m_stagnationIterations;
	LocalSearch *m_opt;
	double lambda;
	double  branching_factor;
	double  trail_0;
	double trail_max;
	double	trail_min;
	bool reset_ant;
	double restart_best_ant_tour_length;

	std::vector<int> best_to_date_tour;  //best to date for even
	std::vector<int> restart_best_ant_path;
	std::vector<int> iteration_best_ant;


	//added
	int m_neg;
	std::vector<std::vector<double> > nn_list;
	void calculateNearestNeigbhor(int NUMOFANTS);
	void init_pheromone_trails(double initial_trail);

	void global_update_pheromone(std::vector<int> &AntTour);

	std::vector<int> m_iteration_best_ant;

 	virtual void updateBestSoFarTour();
	

	int m_bestdistance;
	int m_bestItterdistance;
	void choose_closest_next(std::vector<bool> &antsvisted , std::vector<int> &nntour);
	double node_branching(double l);

public:
	void localSearch();
	void initPheromones();
	void check_pheromone_trail_limits();
	void updatePheromones();
	MinMaxAntSystem(int noNodes, int noAnts, int maxStagnationIterations , std::vector<std::vector<int> >  *Vdistance);
	~MinMaxAntSystem(void);
};

