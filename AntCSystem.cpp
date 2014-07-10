#include "StdAfx.h"
#include "AntCSystem.h"

#include <iostream>


AntCSystem::AntCSystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , MatrixArrayTypeInt *matrix)
					
					: AntSystem(num_cities , ants , *matrix)
{
	  
	this->m_beta = 3;
	this->m_alpha =1;
	this->xi = 0.1;
    this->tau0 = 0;
    this->q0 =0.9;
	this->m_rho = 0.5;
	this->m_r0 = 0.2;
 

}


AntCSystem::~AntCSystem(void)
{
//	delete m_pOpt;
}
void AntCSystem::Evaporation()
{
	for (int i=1; i < m_noNodes; i++) 
	{
		int to = this->m_bestSoFarTour[i];
        int from = this->m_bestSoFarTour[i-1];	
		
		(*m_Vpheromone)[from][to] *= (1.0-m_r0);
		(*m_Vpheromone)[to][from] *= (1.0-m_r0);
		
	}
	
}
/*
In ACS only the globally best ant (i.e., the ant which constructed the shortest tour from the
beginning of the trial) is allowed to deposit pheromone. This choice, together with the use of
the pseudo-random-proportional rule, is intended to make the search more directed: Ants
search in a neighborhood of the best tour found up to the current iteration of the algorithm.
Global updating is performed after all ants have completed their tours. The pheromone level
is updated by applying the global updating rule of Eq. (4)
*/

 void AntCSystem::globalPheromoneDeposit()
{
	//this->get
	double dist = this->getBestSoFarTourLength();// getBestTourLength();
  	double d_tau =  1.0 / m_bestSoFarTourlength;

	for(int i = 0; i < m_noNodes-1; i++)
	{
			int j = m_bestSoFarTour[i];
			int h = m_bestSoFarTour[i+1];

			//(1. - 0.1) * pheromone[h][j] + 0.1 * trail_0
			double value = (*m_Vpheromone)[j][h];
			double updateValue = (1.0 - m_r0) * value + m_r0 * d_tau;
			if(value > updateValue)
				std::cout << "error";

			(*m_Vpheromone)[j][h]= updateValue;
			(*m_Vpheromone)[h][j]= updateValue;
		}
}

void AntCSystem::updateBestSoFarTour()
{
	AntSystem::updateBestSoFarTour();
}
//no need
void AntCSystem::calculateHeuristicMatrix()
{
      for(int i = 0; i < m_noNodes; i++)
          for(int j = 0; j < m_noNodes; j++)
		  {
                double xx = (*m_Vpheromone)[i][j];
				double yy = (*m_Vdistance)[i][j];
				yy = (1.0 / ((double) yy + 0.1));
				(*m_Vchoice)[i][j] = pow( xx , m_alpha) * pow(yy,m_beta);
		  }	
	
	
	
	//AntSystem::calculateHeuristicMatrix();  //same
}
void AntCSystem::updatePheromones()
{
	 globalPheromoneDeposit();
	 Evaporation();
   
	 calculateHeuristicMatrix();
}


void AntCSystem::localPheromoneUpdate(int ant, int step)
{
	int j = this->m_Ants[ant].getTour(step);
    int h = m_Ants[ant].getTour(step-1);

    double currentValue =  (*m_Vpheromone)[h][j];
	double updateValue =  (1.0 - xi)*currentValue+xi*tau0;
    (*m_Vpheromone)[h][j] = updateValue;
    (*m_Vpheromone)[j][h] = updateValue;

  	double dist = (*m_Vdistance)[h][j];
  	double niu = (1.0 / ((double) dist + 0.1)); //stop division by0
    (*m_Vchoice)[h][j] = pow( (*m_Vpheromone)[h][j],m_alpha)*pow(niu,m_beta);
    (*m_Vchoice)[j][h] = (*m_Vchoice)[h][j];
}   
void AntCSystem::initPheromones()
{
	//
	m_bestSoFarTour  = m_opt->greedyTour(this->m_noNodes);
	double lenght = this->calculateTourLength(m_bestSoFarTour);

	this->m_bestSoFarTourlength = (int)lenght;
	tau0= 1.0 / (((double)lenght)*(double)this->m_noAnts);
	for (int i=0; i<m_noNodes; i++)	
	{
		for (int j=0; j<m_noNodes; j++)
			(*m_Vpheromone)[i][j] = tau0;
	}

	for(int i=0;i<this->m_noNodes;i++)
		(*m_Vpheromone)[i][i] =0;


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

void AntCSystem::decisionRule(int k, int step)
{
     //roulette wheel selection  
	int c = m_Ants[k].getTour(step-1); 
 
		for (int i=0; i<m_noNodes; i++) 
		{
			t_prob[i] = 0.0;
			strength[i] = 0;	
			if (m_Ants[k].isVisited(i) == false && i!=c)	
			{
				t_prob[i]= (*m_Vchoice)[c][i];
			}
		}
	//	strength[0] = 0;	
		for (int z =0; z < m_noNodes; z++)
			strength[z+1] = t_prob[z] + strength[z];
				
		double x = fRand(0,  strength[m_noNodes]);
		
		int j = 0;
		while (!((strength[j] <= x) && (x <= strength[j+1])))
		  j++;
		
		int randomDecision =j;
		double maxHeuristic = -1;
        int maxHeuristicIdx = -1;
        for(j = 0; j < m_noNodes; j++)
		{
			 double choice = (*m_Vchoice)[c][j];
			 if(maxHeuristic < choice && !(m_Ants[k].isVisited(j)) && c!=j)
			 {
				maxHeuristic = choice;
				maxHeuristicIdx = j;
			}
        }
		x = fRand(0, 1);	
        if(x < q0)  //if exploitation
		{
            m_Ants[k].setTour(step, maxHeuristicIdx);
            m_Ants[k].setVisited(maxHeuristicIdx);
        }
        else
		{
            m_Ants[k].setTour(step, randomDecision);
            m_Ants[k].setVisited(randomDecision);
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


void AntCSystem::constructSolutions()
{
	//clear ants memory and shuffle start pos
	initAnts();
	 
	//place ants in ramdom citys for starting
	for (int k = 0; k < m_noAnts; k++ )
	{
		m_Ants[k].setTour(0,m_Randomnseq[k]);
		m_Ants[k].setVisited(m_Randomnseq[k]);

	}

	for(int step = 1 ; step < m_noNodes; step++)
	{
		for(int k = 0; k < m_noAnts; k++)
		{
			decisionRule(k,step);
			localPheromoneUpdate(k,step);
		
		}
	}
     
	for(int k = 0; k < m_noAnts; k++)
	{
		//int tourstart=m_Ants[k].getTour(0);
		//m_Ants[k].setTour(m_noNodes,tourstart);
	//	localPheromoneUpdate(k,m_noNodes);
		const std::vector<int>& tourvector = m_Ants[k].getCompleteTour();
		m_Ants[k].settourLength((int)calculateTourLength(tourvector));
	}
	updateBestSoFarTour();
}
