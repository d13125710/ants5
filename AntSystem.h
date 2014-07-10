#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>

#include <algorithm>    // std::random_shuffle


#include "Ant.h"
#include "LocalSearch.h"

typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;
typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class AntSystem
{
		bool m_doLocalSearch , m_doOpt2 , m_doOpt3 , m_doGreed;
protected:

	MatrixArrayTypeDouble *m_Vpheromone;
	MatrixArrayTypeDouble *m_Vchoice;
	MatrixArrayTypeInt  *m_Vdistance;

	LocalSearch *m_opt;

	std::mt19937   m_generator;
	
	std::vector<Ant> m_Ants;
	std::vector<int> m_bestSoFarTour;
	std::vector<int> m_Randomnseq;
	std::vector<int> m_bestIrrTour;

	std::vector<double> t_prob;
	std::vector<double> strength;

			
    int m_noNodes;
    double m_alpha;
    double m_beta;
    double m_rho;
    int m_iterations;
	int m_noAnts;
	int m_bestSoFarTourlength;

	void Create(int noNodes, int noAnts);

public:

	virtual void initPheromones();
	void initAnts();
	void initMatrix(MatrixArrayTypeInt &distanceMatrix);
	AntSystem(int noNodes, int noAnts);
	AntSystem(int noNodes, int noAnts, MatrixArrayTypeInt &distanceMatrix);

	std::vector<int> &getBestTourIrreration();
	double calculateTourLength(const std::vector<int> &currentTour);
	
	virtual void Evaporation();
	void adjust_pheremone(int antPos);
	void decisionRule(int k, int step);

	//virtual
	virtual void constructSolutions();
	virtual void updateBestSoFarTour();
	virtual void updatePheromones();
	virtual void calculateHeuristicMatrix();
	
	virtual ~AntSystem(void);

	int getIteration() const
	{
		return m_iterations;
	}

	int getItterationBestTourLength()
	{
		int best= std::numeric_limits<int>::max();
		for( int i=0; i<  m_noAnts ; i++)
		{
			int length = m_Ants[i].getAntTourLength();
			if( length < best)
				best =length;

		}
		
		//std::vector<int> tour=getBestTourIrreration(); //const ref binds to temp objects
		return best;
	}
	double getBestSoFarTourLength()
	{
		return calculateTourLength(m_bestSoFarTour);
	}
	void incIteration()
	{
		m_iterations++;
	}
	int  getIteration()
	{
		return m_iterations;
	}

	const std::vector<int> &getBestSoFarTour()
	{
		return m_bestSoFarTour;
	}
	double getBestTourlengthSofar()
	{
		return m_bestSoFarTourlength;
	}
	void setPheromones(std::vector<std::vector<double> > &MatrixArrayType)
	{
		this->m_Vpheromone = &MatrixArrayType;
	}
	void setDistance(std::vector<std::vector<int> > &MatrixArrayType)
	{
		this->m_Vdistance = &MatrixArrayType;
	}
	const Ant& getAnt(int pos)
	{
      return m_Ants[pos];
    }
	double getPheromone(int i, int j) const
	{
        return (*m_Vpheromone)[i][j];
    }
	double getchoice(int i, int j) const
	{
        return (*m_Vchoice)[i][j];
    }
    int getDistance(int i, int j) const
	{
        return (*m_Vdistance)[i][j];
    }

	void localSearch()
{

        if(m_doOpt2){
            for(int i = 0; i < m_noAnts; i++){
			   m_opt->opt2(m_Ants[i].getCompleteTour());
            }
        }
        if(m_doOpt3){
            for(int i = 0; i < m_noAnts; i++){
                m_opt->opt3(m_Ants[i].getCompleteTour());
            }
        }
		if(m_doGreed){
            for(int i = 0; i < m_noAnts; i++){
                m_opt->greedyAntTour(m_Ants[i].getCompleteTour());
            }
        }


        updateBestSoFarTour();
    
    }


	
};

