#include "StdAfx.h"
#include "CMinMaxAntSystem.h"

#include <iostream>

CMinMaxAntSystem::CMinMaxAntSystem(unsigned int noNodes, unsigned int noAnts, unsigned int maxStagnationIterations , std::vector<std::vector<int> >  *Vdistance) 
	: CAntSystem(noNodes , noAnts , Vdistance)
{
	 m_stagnationIterations = 0;
     m_maxStagnationIterations = maxStagnationIterations;
     m_doOpt2 = false;
     m_doOpt3 = true;
	 //m_pLocalSearch = new CLocalSearch(Vdistance);
//	 m_bestDistanceLength = std::numeric_limits<int>::max();
	 m_bestIterationLength = std::numeric_limits<int>::max();
	 m_restartBestAntTourLength =  std::numeric_limits<int>::max();
	 lambda = 0.05;
	 m_resetAnt =false;
//	 m_bestToDatePath.resize(noNodes);  //best t o date for even
	 m_restartAntBestPath.resize(noNodes); ;
	 m_bestIterationAntPath.resize(noNodes); ;

}
void CMinMaxAntSystem::chooseClosestNext(std::vector<bool> &antsvisted , std::vector<int> &nntour)
{
	unsigned int city, current_city, next_city;
	int min_distance;
	next_city =this->m_noNodes;

	for(unsigned int phase = 1; phase < m_noNodes; phase++)
	{
		current_city = nntour[phase - 1];
		min_distance = std::numeric_limits<int>::max();; /* Search shortest edge */
		
			for (city = 0; city < m_noNodes; city++)
			{
				if (antsvisted[city])
				; /* city already visited */
				else
				{
				
					if ((*m_distanceMatrix)[current_city][city] < min_distance)
					{
						next_city = city;
						min_distance = (*m_distanceMatrix)[current_city][city];
					}
				}
			}
	
	
		nntour[phase] = next_city;
		antsvisted[next_city] = true;
	}

}

 double CMinMaxAntSystem::nodeBranching(double l)
 {
	return this->m_newPheromoneMatrix->average_lambda_branching_factor(l);
}

void CMinMaxAntSystem::calculateNearestNeigbhor(unsigned int NUMOFANTS)
{
	unsigned int i = 0;
	int *pHelpArray = new int[m_noNodes];
	double *pdistanceArray = new double[m_noNodes ];
	*(pdistanceArray) = std::numeric_limits<double>::max();

	m_nnList.resize(m_noNodes);
		for (unsigned int i = 0; i < m_noNodes; ++i)
			m_nnList[i].resize(NUMOFANTS);

	
	for (unsigned int node = 0; node < m_noNodes; node++) 
	{ 

			for (i = 0; i < m_noNodes; i++) 
			{
				*(pdistanceArray+i)=(*m_distanceMatrix)[node][i];
				*(pHelpArray+i)=  i;
			}
			double max = std::numeric_limits<double>::max() - 1; 
			*(pdistanceArray+node) =  std::numeric_limits<double>::max();  // set to a large value .. 
			this->m_pLocalSearch->sort2(pdistanceArray, pHelpArray, 0, static_cast<unsigned int>(m_noNodes - 1));
			for (i = 0; i < NUMOFANTS; i++) 
			{
				m_nnList[node][i] = *(pHelpArray+i);
			}
	}
	delete [] pHelpArray;
	delete [] pdistanceArray;
}

void CMinMaxAntSystem::updateBestSoFarPath()
{
		CAntSystem::updateBestSoFarPath();
		//best in itteration
		if (m_bestSoFarPathlength < m_restartBestAntTourLength)
		{
				m_restartBestAntTourLength =m_bestSoFarPathlength;
				m_restartAntBestPath = getBestTourIrreration();
				irreationSinceLastBest =0;
		}
		else 
		{
			m_bestIterationLength = m_bestSoFarPathlength;
			m_bestIterationAntPath = this->getBestTourIrreration();
			m_branchingFactor = this->nodeBranching(lambda);
			double p_x = exp(log(0.05) / m_noNodes);
			trail_min = 1. * (1. - p_x) / (p_x * (double) ((this->m_noNodes + 1) / 2));
			trail_max = 1. / ((m_rho) * m_bestIterationLength);
			trail_0 = trail_max;
			trail_min = trail_max * trail_min;
			irreationSinceLastBest ++;
	
        }
 }
void CMinMaxAntSystem::initPheromones()
{
		std::vector<bool> visited(m_noNodes);
		std::vector<int> nntour(m_noNodes);
		calculateNearestNeigbhor(20);
		//calculate min max values inital
		int phase = 0;
		int rnd= (rand()%(visited.size()-1))+1;
		nntour[0] =rnd;
		visited[rnd] = true;
		chooseClosestNext(visited, nntour);
		double distance = this->calculatePathLength( nntour);
		trail_max = 1. / ((m_rho) * distance);
		trail_min = trail_max / (2. * m_noNodes);
		initPheromoneTrails(trail_max);
}	
	
//There are variants in the selection of the ants allowed to update pheromones: the best to date ant, 
//or the best for current iteration, or the best after latest reset ant, or the best to date ant for even iterations, 
//and the best for iteration for odd iterations.

void CMinMaxAntSystem::initPheromoneTrails(double initialValue) const
 {
	unsigned int i, j;

		for (i = 0; i < m_noNodes; i++)
		{
			for (j = 0; j <= i; j++) 
			{
				m_newPheromoneMatrix->set(i , j , initialValue);
			//	m_newPheromoneMatrix->set(j , i , initialValue);
				
			//	(*m_pheromoneMatrix)[i][j] = initialValue;
			//	(*m_pheromoneMatrix)[j][i] = initialValue;
			}
		}

		for(i = 0; i < m_noNodes; i++)
			m_newPheromoneMatrix->set(i , i , 0);
}
void CMinMaxAntSystem::globalUpdatePheromone(const std::vector<int> &AntTour)
{
	int Q=1;


	double tour_length = calculatePathLength (AntTour); 
	

	for(int city = 1; city < AntTour.size(); city++)
	 {
	      int from = AntTour[city-1];
		  int to = AntTour[city];
 
		  // eq 14.2 / 14.3
		 // (*m_pheromoneMatrix)[from][to] += (Q/tour_length) ;
		  //(*m_pheromoneMatrix)[to][from] += (Q/tour_length) ; 
		  double tpadd = (Q/tour_length);
		  m_newPheromoneMatrix->add(from , to , tpadd);
	 	 //  m_newPheromoneMatrix->add(to , from , tpadd );
		  //m_newPheromoneMatrix->add(to , from , (Q/tour_length) );
	
  }


	

}
void CMinMaxAntSystem::checkPheromoneLimits()
{
	unsigned int  i, j;

		for (i = 0; i < m_noNodes ; i++) 
		{
			for (j = 0; j < i; j++) 
			{
				double value =  m_newPheromoneMatrix->get(i,j);
				 
				if ( value < trail_min)
				{
					m_newPheromoneMatrix->set(i,j , trail_min);
				//	m_newPheromoneMatrix->set(j,i, trail_min);
					//(*m_pheromoneMatrix)[i][j] = trail_min;
					//(*m_pheromoneMatrix)[j][i] = trail_min;
				} 
				else if ( value > trail_max )
				{
					m_newPheromoneMatrix->set(i,j , trail_max);
				//	m_newPheromoneMatrix->set(j,i , trail_max);
				//	(*m_pheromoneMatrix)[i][j] = trail_max;
					//(*m_pheromoneMatrix)[j][i] = trail_max;
				}


			/*	if ( (*m_pheromoneMatrix)[i][j] < trail_min)
				{
					 (*m_pheromoneMatrix)[i][j] = trail_min;
					 (*m_pheromoneMatrix)[j][i] = trail_min;
				} 
				else if ( (*m_pheromoneMatrix)[i][j] > trail_max)
				{
					 (*m_pheromoneMatrix)[i][j] = trail_max;
					 (*m_pheromoneMatrix)[j][i] = trail_max;
				}
				*/
			}
		}
}


void CMinMaxAntSystem::updatePheromones()
{
	//end of ants walk
	evaporateAllPhero ();	

	
	if (this->m_iterations % 25 == 0) //even
	{
	  
		//std::cout << "globalUpdatePheromone best so far " << calculatePathLength(m_bestSoFarPath) << ", iteration: "  << m_iterations << ", b_fac " 
		//	<< m_branchingFactor << std::endl;

		globalUpdatePheromone(this->m_bestSoFarPath);  //best to date for even
	} 
	else //odd
	{
	    if(m_resetAnt)
		{
			globalUpdatePheromone(m_restartAntBestPath);
			m_resetAnt = false;
		}
	    else
		{
			std::vector<int> &best = getBestTourIrreration();
			globalUpdatePheromone(best);
	//		std::cout << "globalUpdatePheromone best for ittteration " << calculatePathLength(best) << ", iteration: "  << m_iterations << ", b_fac " 
	//			<< m_branchingFactor << std::endl;

		}
	}
	

			
	if ((m_iterations % 100) == 0 && irreationSinceLastBest > 150) 
	{
		

		
		//create nn_list ..
		m_branchingFactor = nodeBranching(lambda);
		double  branch_fac = 1.00001;
		std::cout << "best so far " << calculatePathLength(m_bestSoFarPath) << ", iteration: "  << m_iterations << ", b_fac " << m_branchingFactor << "avg" << std::endl;

		if ( (m_branchingFactor < branch_fac)) 
		{
			std::cout << "INIT TRAILS!!!\n";
			//m_restartBestAntTourLength
			m_restartBestAntTourLength =  std::numeric_limits<int>::max();
			initPheromoneTrails(trail_max);
			m_resetAnt = true;
			irreationSinceLastBest =0;

		}

	}

	checkPheromoneLimits();

	this->calculateHeuristicMatrix();



    
}
/*
void CMinMaxAntSystem::localSearch()
{
        if(m_doOpt2){
            for(unsigned int i = 0; i < m_noAnts; i++){
			   m_pLocalSearch->opt2(m_Ants[i].getAntsCityTour());
            }
        }
        if(m_doOpt3){
            for(unsigned int i = 0; i < m_noAnts; i++){
                m_pLocalSearch->opt3(m_Ants[i].getAntsCityTour());
            }
        }

        updateBestSoFarTour();
        if(m_stagnationIterations > 0)
           m_stagnationIterations--;
    }

*/

CMinMaxAntSystem::~CMinMaxAntSystem(void)
{
	//if(m_pLocalSearch)
	//	delete m_pLocalSearch;
}
