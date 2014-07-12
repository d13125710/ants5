#include "StdAfx.h"
#include "CAntColonySystem.h"

#include <iostream>


CAntColonySystem::CAntColonySystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , MatrixArrayTypeInt *matrix)
					
					: CAntSystem(num_cities , ants , matrix)
{
	  
	this->m_beta = 3;
	this->m_alpha =1;
	this->xi = 0.1;
    this->tau0 = 0;
    this->q0 =0.9;
	this->m_rho = 0.5;
	this->m_r0 = 0.1;
 

}

CAntColonySystem::~CAntColonySystem(void)
{
//	delete m_pOpt;
}
void CAntColonySystem::evaporateAllPhero()
{
	for (unsigned int i=1; i < m_noNodes; i++) 
	{
		int to = this->m_bestSoFarPath[i];
        int from = this->m_bestSoFarPath[i-1];	
		
		(*m_pheromoneMatrix)[from][to] *= (1.0-m_r0);
		(*m_pheromoneMatrix)[to][from] *= (1.0-m_r0);
		
	}
	//end phero try
//	(*m_Vpheromone)[m_noNodes-1][0] *=  1.0-m_r0;

	
}
/*
In ACS only the globally best ant (i.e., the ant which constructed the shortest tour from the
beginning of the trial) is allowed to deposit pheromone. This choice, together with the use of
the pseudo-random-proportional rule, is intended to make the search more directed: Ants
search in a neighborhood of the best tour found up to the current iteration of the algorithm.
Global updating is performed after all ants have completed their tours. The pheromone level
is updated by applying the global updating rule of Eq. (4)
*/

 void CAntColonySystem::globalPheromoneDeposit()
{
	//this->get
	double dist = this->getBestSoFarPathLength();// getBestTourLength();
  	double d_tau =  1.0 / m_bestSoFarPathlength;

	for(unsigned int i = 0; i < m_noNodes-1; i++)
	{
			int j = m_bestSoFarPath[i];
			int h = m_bestSoFarPath[i+1];

			//(1. - 0.1) * pheromone[h][j] + 0.1 * trail_0
			double value = (*m_pheromoneMatrix)[j][h];
			double updateValue = (1.0 - m_r0) * value + m_r0 * d_tau;
			if(value > updateValue)
				std::cout << "error";

			(*m_pheromoneMatrix)[j][h]= updateValue;
			(*m_pheromoneMatrix)[h][j]= updateValue;
		}
	//update the end pheromone try
	//(*m_Vpheromone)[m_noNodes-1][0]=  (1.0 - m_r0) * (*m_Vpheromone)[m_noNodes-1][0] + m_r0 * d_tau;

}

void CAntColonySystem::updateBestSoFarTour()
{
	CAntSystem::updateBestSoFarPath();
}
//no need
void CAntColonySystem::calculateHeuristicMatrix()
{
      for(unsigned int i = 0; i < m_noNodes; i++)
          for(unsigned int j = 0; j < m_noNodes; j++)
		  {
                double xx = (*m_pheromoneMatrix)[i][j];
				double yy = (*m_distanceMatrix)[i][j];
				yy = (1.0 / ((double) yy + 0.1));
				(*m_heuristicMatrix)[i][j] = pow( xx , m_alpha) * pow(yy,m_beta);
		  }	
	
	
	
	//AntSystem::calculateHeuristicMatrix();  //same
}
void CAntColonySystem::updatePheromones()
{
	evaporateAllPhero();
	globalPheromoneDeposit();
	calculateHeuristicMatrix();
}


void CAntColonySystem::localPheromoneUpdate(int ant, int step)
{
	int j = this->m_Ants[ant].getCity(step);
    int h = m_Ants[ant].getCity(step-1);

    double currentValue =  (*m_pheromoneMatrix)[h][j];
	double updateValue =  (1.0 - xi)*currentValue+xi*tau0;
    (*m_pheromoneMatrix)[h][j] = updateValue;
    (*m_pheromoneMatrix)[j][h] = updateValue;

  	double dist = (*m_distanceMatrix)[h][j];
  	double niu = (1.0 / ((double) dist + 0.1)); //stop division by0
    (*m_heuristicMatrix)[h][j] = pow( (*m_pheromoneMatrix)[h][j],m_alpha)*pow(niu,m_beta);
    (*m_heuristicMatrix)[j][h] = (*m_heuristicMatrix)[h][j];
}   
void CAntColonySystem::initPheromones()
{
	//
	m_bestSoFarPath  = m_pLocalSearch->greedyPath(this->m_noNodes);
	double lenght = this->calculatePathLength(m_bestSoFarPath);

	this->m_bestSoFarPathlength = (int)lenght;
	tau0= 1.0 / (((double)lenght)*(double)this->m_noAnts);
	for (unsigned int i=0; i<m_noNodes; i++)	
	{
		for (unsigned int j=0; j<m_noNodes; j++)
			(*m_pheromoneMatrix)[i][j] = tau0;
	}

	for(unsigned int i=0;i<this->m_noNodes;i++)
		(*m_pheromoneMatrix)[i][i] =0;


}
/*
The procedure works as follows: first, the current city c of ant k is determined
(line 1). The probabilistic choice of the next city then works analogously
to the roulette wheel selection procedure of evolutionary computation (Goldberg,
1989): each value choice_info½c½ j of a city j that ant k has not visited yet determines
a slice on a circular roulette wheel, the size of the slice being proportional to the
weight of the associated choice (lines 2–10). Next, the wheel is spun and the city to
which the marker points is chosen as the next city for ant k (lines 11–17). This is
implemented by
1. summing the weight of the various choices in the variable sum_probabilities,
2. drawing a uniformly distributed random number r from the interval ½0; sum_
probabilities,
3. going through the feasible choices until the sum is greater or equal to r.
Finally, the ant is moved to the chosen city, which is marked as visited (lines 18
and 19).
*/

void CAntColonySystem::decisionRule(int k, int step)
{
     //roulette wheel selection  
	int c = m_Ants[k].getCity(step-1); 
 
		for (unsigned int i=0; i<m_noNodes; i++) 
		{
			t_prob[i] = 0.0;
			m_strength[i] = 0; 	
			if (m_Ants[k].isCityVisited(i) == false)// && i!=c)	
			{
				t_prob[i]= (*m_heuristicMatrix)[c][i];
			}
		}
	//	strength[0] = 0;	
		for (unsigned int z =0; z < m_noNodes; z++)
			m_strength[z+1] = t_prob[z] + m_strength[z];
				
		double x = fRand(0,  m_strength[m_noNodes]);
		
		unsigned int j = 0;
		while (!((m_strength[j] <= x) && (x <= m_strength[j+1])))
		  j++;
		
		int randomDecision =j;
		double maxHeuristic = -1;
        int maxHeuristicIdx = -1;
        for(j = 0; j < m_noNodes; j++)
		{
			 double choice = (*m_heuristicMatrix)[c][j];
			 if(maxHeuristic < choice && !(m_Ants[k].isCityVisited(j)))// && c!=j)
			 {
				maxHeuristic = choice;
				maxHeuristicIdx = j;
			}
        }
		x = fRand(0, 1);	
        if(x < q0)  //if exploitation
		{
            m_Ants[k].setAntCity(step, maxHeuristicIdx);
            m_Ants[k].setCityVisited(maxHeuristicIdx);
        }
        else
		{
			//	std::cout << maxHeuristicIdx << ",";
            m_Ants[k].setAntCity(step, randomDecision);
            m_Ants[k].setCityVisited(randomDecision);
        }

}
/*
1. First, the ants’ memory must be emptied. This is done in lines 1 of procedure
ConstructSolutions by marking all cities as unvisited, that is, by setting all the entries
of the array ants.visited to false for all the ants.
2. Second, each ant has to be assigned an initial city. One possibility is to assign each
ant a random initial city. This is accomplished in lines 6 to 11 of the procedure. The
function random returns a random number chosen according to a uniform distribution
over the set f1; . . . ; ng.
3. Next, each ant constructs a complete tour. At each construction step (see the procedure
in figure 3.17) the ants apply the AS action choice rule [equation (3.2)]. The
procedure ASDecisionRule implements the action choice rule and takes as parameters
the ant identifier and the current construction step index; this is discussed below in
more detail.
4. Finally, in lines 18 to 21, the ants move back to the initial city and the tour length
of each ant’s tour is computed. Remember that, for the sake of simplicity, in the tour
representation we repeat the identifier of the first city at position n þ 1; this is done in
line 19.

*/


void CAntColonySystem::constructSolutions()
{
	//clear ants memory and shuffle start pos
	initAnts();
	 
	//place ants in ramdom citys for starting
	for (unsigned int k = 0; k < m_noAnts; k++ )
	{
		m_Ants[k].setAntCity(0,m_randomPath[k]);
		m_Ants[k].setCityVisited(m_randomPath[k]);

	}

	for(unsigned int step = 1 ; step < m_noNodes; step++)
	{
		for(unsigned int k = 0; k < m_noAnts; k++)
		{
			decisionRule(k,step);
			localPheromoneUpdate(k,step);
		
		}
	}
     
	for(unsigned int k = 0; k < m_noAnts; k++)
	{
		int tourstart=m_Ants[k].getCity(0);
		m_Ants[k].setAntCity(m_noNodes,tourstart);
		localPheromoneUpdate(k,m_noNodes);
		const std::vector<int>& tourvector = m_Ants[k].getAntsCityTour();
		m_Ants[k].setAntTourLength((int)calculatePathLength(tourvector));
	}
	updateBestSoFarTour();
}
