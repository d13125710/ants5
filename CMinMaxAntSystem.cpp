#include "StdAfx.h"
#include "CMinMaxAntSystem.h"


CMinMaxAntSystem::CMinMaxAntSystem(int noNodes, int noAnts, int maxStagnationIterations , std::vector<std::vector<int> >  *Vdistance) 
	: CAntSystem(noNodes , noAnts , Vdistance)
{
	 m_stagnationIterations = 0;
     m_maxStagnationIterations = maxStagnationIterations;
     m_doOpt2 = false;
     m_doOpt3 = true;
	 m_pLocalSearch = new CLocalSearch(Vdistance);
	 m_bestDistanceLength = std::numeric_limits<int>::max();
	 m_bestIterationLength = std::numeric_limits<int>::max();
	 m_restartBestAntTourLength =  std::numeric_limits<int>::max();
	
	 m_neg = 20;
	  if(m_neg > noNodes)
		 m_neg = noNodes-1;
	 lambda = 0.05;
	 m_resetAnt =false;
	 m_bestToDatePath.resize(noNodes);  //best to date for even
	 m_restartAntBestPath.resize(noNodes); ;
	 m_bestIterationAntPath.resize(noNodes); ;

}
void CMinMaxAntSystem::chooseClosestNext(std::vector<bool> &antsvisted , std::vector<int> &nntour)
{
	int city, current_city, next_city, min_distance;
	next_city =this->m_noNodes;

	for(int phase = 1; phase < m_noNodes; phase++)
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
	int i, m;
	double min, max, cutoff;
	double avg;

	std::vector<double> num_branches(m_noNodes);

	
	for (m = 0; m < m_neg; m++) 
	{
	   
	    min = (*m_pheromoneMatrix)[m][m_nnList[m][1]];
	    max = (*m_pheromoneMatrix)[m][m_nnList[m+1][1]];


	    for (i = 1; i < m_neg; i++) 
		{
			if ((*m_pheromoneMatrix)[m][m_nnList[m][i]] > max)
				max = (*m_pheromoneMatrix)[m][m_nnList[m][i]];
			if ((*m_pheromoneMatrix)[m][m_nnList[m][i]] < min)
				min = (*m_pheromoneMatrix)[m][m_nnList[m][i]];
	    }
	    cutoff = min + l * (max - min);

	    for (i = 1; i < m_neg; i++)
		{
			if ((*m_pheromoneMatrix)[m][m_nnList[m][i]] > cutoff)
				num_branches[m] += 1;
	    }
	}
	avg = 0.;
	for (m = 0; m < m_noNodes; m++) {
	    avg += num_branches[m];
	}
	
	return (avg / (double) (m_noNodes * 2));
}

void CMinMaxAntSystem::calculateNearestNeigbhor(int NUMOFANTS)
{
	int i = 0;
	double *pHelpArray = new double[m_noNodes];
	double *pdistanceArray = new double[m_noNodes ];
	*(pdistanceArray) = std::numeric_limits<double>::max();

	m_nnList.resize(m_noNodes);
		for (int i = 0; i < m_noNodes; ++i)
			m_nnList[i].resize(NUMOFANTS);

	
	for (int node = 0; node < m_noNodes; node++) 
	{ 

			for (i = 0; i < m_noNodes; i++) 
			{
				*(pdistanceArray+i)=(*m_distanceMatrix)[node][i];
				*(pHelpArray+i)=  i;
			}
			double max = std::numeric_limits<double>::max() - 1; 
			*(pdistanceArray+node) =  std::numeric_limits<double>::max();  // set to a large value .. 
			this->m_pLocalSearch->sort2(pdistanceArray, pHelpArray, 0, m_noNodes - 1);
			for (i = 0; i < NUMOFANTS; i++) 
			{
				m_nnList[node][i] = *(pHelpArray+i);
			}
	}
	delete [] pHelpArray;
	delete [] pdistanceArray;
}

void CMinMaxAntSystem::updateBestSoFarTour()
{
	
	    m_bestSoFarPathlength= (int) calculatePathLength(m_bestSoFarPath);

		
		if(this->m_iterations == 0) //first init
		{
			m_bestSoFarPath = getBestTourIrreration();
			m_bestSoFarPathlength= (int) calculatePathLength(m_bestSoFarPath);
			m_stagnationIterations =0;
		}

		int d = (int) getItterationBestPathLength() ;

		//best in itteration
		if (d < m_restartBestAntTourLength)
		{
				m_restartBestAntTourLength = d;
				m_restartAntBestPath = getBestTourIrreration();
				m_stagnationIterations =0;
		}

		if(d < m_bestSoFarPathlength)
		{
            m_bestDistanceLength = d;

			m_bestToDatePath= getBestTourIrreration();  //best to date for even
			m_bestSoFarPath = getBestTourIrreration();
			m_stagnationIterations = 0;
        }
		else 
		{
            //copy the tour
			//iteration_best_ant
			m_bestIterationLength = d;
			m_iterationBestAntPath = getBestTourIrreration();
			m_branchingFactor = nodeBranching(lambda);
			double p_x = exp(log(0.05) / m_noNodes);
			trail_min = 1. * (1. - p_x) / (p_x * (double) ((this->m_noAnts + 1) / 2));
			trail_max = 1. / ((m_rho) * m_bestIterationLength);
			trail_0 = trail_max;
			trail_min = trail_max * trail_min;
			m_stagnationIterations ++;
	
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
	int i, j;

		for (i = 0; i < m_noNodes; i++)
		{
			for (j = 0; j <= i; j++) 
			{
				(*m_pheromoneMatrix)[i][j] = initialValue;
				(*m_pheromoneMatrix)[j][i] = initialValue;
			}
		}
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
		  (*m_pheromoneMatrix)[from][to] += (Q/tour_length) ;
		  (*m_pheromoneMatrix)[to][from] += (Q/tour_length) ; 
	
  }

}
void CMinMaxAntSystem::checkPheromoneLimits()
{
	int i, j;

		for (i = 0; i < m_noNodes ; i++) 
		{
			for (j = 0; j <= i; j++) 
			{
				if ( (*m_pheromoneMatrix)[i][j] < trail_min)
				{
					 (*m_pheromoneMatrix)[i][j] = trail_min;
					 (*m_pheromoneMatrix)[j][i] = trail_min;
				} 
				else if ( (*m_pheromoneMatrix)[i][j] > trail_max)
				{
					 (*m_pheromoneMatrix)[i][j] = trail_max;
					 (*m_pheromoneMatrix)[j][i] = trail_max;
				}
			}
		}
}


void CMinMaxAntSystem::updatePheromones()
{
	evaporateAllPhero ();	
	if (this->m_iterations % 50 == 0) //even
	{
	  
		globalUpdatePheromone(m_bestToDatePath);  //best to date for even
	} 
	else //odd
	{
	    if(m_resetAnt)
			globalUpdatePheromone(m_restartAntBestPath);
	    else
			globalUpdatePheromone(getBestTourIrreration());
	}
	
    //end of ants walk
	
	checkPheromoneLimits();
			
	if ((m_iterations % 100) == 0 && m_stagnationIterations > 20) 
	{
		//create nn_list ..
		m_branchingFactor = nodeBranching(lambda);
		double  branch_fac = 1.00001;
		if ( (m_branchingFactor < branch_fac)	&& (m_iterations  > 150)) 
		{
			//cout << "INIT TRAILS!!!\n";
			m_restartBestAntTourLength =  std::numeric_limits<int>::max();
			initPheromoneTrails(trail_max);
			m_resetAnt = true;
			m_stagnationIterations =0;

		}

	}
	else
	{
		m_resetAnt = false;
	}

	this->calculateHeuristicMatrix();



    
}

void CMinMaxAntSystem::localSearch()
{
        if(m_doOpt2){
            for(int i = 0; i < m_noAnts; i++){
			   m_pLocalSearch->opt2(m_Ants[i].getAntsCityTour());
            }
        }
        if(m_doOpt3){
            for(int i = 0; i < m_noAnts; i++){
                m_pLocalSearch->opt3(m_Ants[i].getAntsCityTour());
            }
        }

        updateBestSoFarTour();
        if(m_stagnationIterations > 0)
           m_stagnationIterations--;
    }



CMinMaxAntSystem::~CMinMaxAntSystem(void)
{
	if(m_pLocalSearch)
		delete m_pLocalSearch;
}
