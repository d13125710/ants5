#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>

#include <algorithm>    // std::random_shuffle


#include "CAnt.h"
#include "CLocalSearch.h"

typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;
typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class CAntSystem
{
	bool m_bLocalSearch , m_LocalSearchOpt2 , m_LocalSearchOpt3 , m_LocalSearchGreed;

	

protected:

	MatrixArrayTypeDouble *m_pheromoneMatrix;
	MatrixArrayTypeDouble *m_heuristicMatrix;
	MatrixArrayTypeInt  *m_distanceMatrix;
	long irreationSinceLastBest;
	CLocalSearch *m_pLocalSearch;

//	std::mt19937   m_randomGenerator;
	
	std::vector<CAnt> m_Ants;
	std::vector<int> m_bestSoFarPath;
	std::vector<int> m_randomPath;
	std::vector<int> m_bestIrrPath;
	std::vector<double> t_prob;
	std::vector<double> m_strength;

			
    int m_noNodes;
    double m_alpha;
    double m_beta;
    double m_rho;
    int m_iterations;
	int m_noAnts;
	int m_bestSoFarPathlength;
	int BestAntIndex;

	void Create(int noNodes, int noAnts);
	double fRand(double a, double b)
	{
		static bool seeded = false;
		if(!seeded) {
			srand(time(NULL));
			seeded = true;
		}
		long _rnd = rand();
		if (_rnd == 0)
			_rnd++;
		return b * (double)_rnd / RAND_MAX;	
	}


public:

	virtual void initPheromones();
	void initAnts();
	void initMatrix(MatrixArrayTypeInt &distanceMatrix);
	CAntSystem(int noNodes, int noAnts);
	CAntSystem(int noNodes, int noAnts, MatrixArrayTypeInt *distanceMatrix);

	std::vector<int> &getBestTourIrreration();
	double calculatePathLength(const std::vector<int> &currentPath) const;
	
	virtual void evaporateAllPhero();
	void adjustPheremone(int antPos);
	void decisionRule(int k, int step);

	//virtual
	virtual void constructSolutions();
	virtual void updateBestSoFarPath();
	virtual void updatePheromones();
	virtual void calculateHeuristicMatrix();
	
	virtual ~CAntSystem(void);

	int getIteration() const
	{
		return m_iterations;
	}

	int getItterationBestPathLength() const
	{
		int best= std::numeric_limits<int>::max();
		for( int i=0; i<  m_noAnts ; i++)
		{
			int length = m_Ants[i].getAntTourLength();
			if( length < best)
				best =length;

		}
		return best;
	}
	double getBestSoFarPathLength() const
	{
		return calculatePathLength(m_bestSoFarPath);
	}
	void incIteration() 
	{
		m_iterations++;
	}
	const std::vector<int> &getBestSoFarPath()
	{
		return m_bestSoFarPath;
	}
	double getBestPathLengthSofar()
	{
		return m_bestSoFarPathlength;
	}
	void setPheromonesMatrix(std::vector<std::vector<double> > &MatrixArrayType)
	{
		this->m_pheromoneMatrix = &MatrixArrayType;
	}
	void setDistance(std::vector<std::vector<int> > &MatrixArrayType)
	{
		this->m_distanceMatrix = &MatrixArrayType;
	}
	const CAnt& getAnt(int pos)
	{
      return m_Ants[pos];
    }
	const CAnt &GetBestAnt()
	{
		return getAnt(BestAntIndex);
	}
	double getPheromone(int i, int j) const
	{
        return (*m_pheromoneMatrix)[i][j];
    }
	double getHeuristicAt(int i, int j) const
	{
        return (*m_heuristicMatrix)[i][j];
    }
    int getDistance(int i, int j) const
	{
        return (*m_distanceMatrix)[i][j];
    }

	void localSearch()
{

        if(m_LocalSearchOpt2){
            for(int i = 0; i < m_noAnts; i++){
			   m_pLocalSearch->opt2(m_Ants[i].getAntsCityTour());
            }
        }
        if(m_LocalSearchOpt3){
            for(int i = 0; i < m_noAnts; i++){
                m_pLocalSearch->opt3(m_Ants[i].getAntsCityTour());
            }
        }
		if(m_LocalSearchGreed){
            for(int i = 0; i < m_noAnts; i++){
                m_pLocalSearch->greedyAntPath(m_Ants[i].getAntsCityTour());
            }
        }


        updateBestSoFarPath();
    
    }


	
};

