#include "StdAfx.h"
#include "CAntColonySystem.h"

#include <iostream>




CAntColonySystem::CAntColonySystem(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , MatrixArrayTypeInt *matrix)
					
					: CAntSystem(num_cities , ants , matrix)
{
	  
	this->m_beta = 2;
	this->m_alpha =1;
	this->xi = 0.1;
    this->tau0 = 0;
    this->q0 =0.90;
	this->m_rho = 0.5;
	this->m_r0 = 0.3;
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
		double value = m_newPheromoneMatrix->get(from,to);
		
		m_newPheromoneMatrix->set(from,to, (1. - m_rho) * value + m_rho * this->m_bestSoFarPathlength);

	//	m_newPheromoneMatrix->multipy(from , to ,(1-m_rho));
	////	pheromone[idx1][idx2]*=(1-ro);
//this->m_newPheromoneMatrix->evaporate(from , to);
		//this->m_newPheromoneMatrix->evaporate(to , from);
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

 void CAntColonySystem::globalPheromoneDeposit()
{
	//this->get
	double dist = m_bestSoFarPathlength;
  	double d_tau =  1.0 / dist;

	for(unsigned int i = 1; i < m_noNodes; i++)
	{ 
		int to = m_bestSoFarPath[i]; //h
		int from = m_bestSoFarPath[i-1]; //j
		double value = this->m_newPheromoneMatrix->get(from , to);
		double updateValue = (1.0 - m_r0) * value + m_r0 * d_tau;
		//if(value > updateValue) 
		//	std::cout << "error";
		this->m_newPheromoneMatrix->set( from , to , updateValue);


		double dist = (*m_distanceMatrix)[to][from];
		(*m_heuristicMatrix)[from][to] = pow( updateValue,m_alpha) * pow(HEURISTIC(from,to),m_beta);
		(*m_heuristicMatrix)[to][from] = (*m_heuristicMatrix)[from][to];
	}
	
}

void CAntColonySystem::updateBestSoFarPath()
{
	//find best ant
	CAnt *pAnt=0;
	int best_distance=999999;
	for(unsigned int i = 0; i < this->m_Ants.size(); i++)
	{
		int currentL = (int) calculatePathLength(m_Ants[i].getAntsCityTour());
		if(currentL < best_distance)
		{
			best_distance = currentL;
			pAnt = &m_Ants[i];

		}
	}
	//pAmt is the best ant;
	
	if(best_distance < m_bestSoFarPathlength )
	{

		for(unsigned int j = 0; j < m_noNodes; j++)
		{
			m_bestSoFarPath[j] = pAnt->getCity(j);

		}

		m_bestSoFarPathlength = best_distance;

	}


	//CAntSystem::updateBestSoFarPath();
}
//no need
void CAntColonySystem::calculateHeuristicMatrix()
{
   // CAntSystem::calculateHeuristicMatrix();
	

	for(unsigned int i = 0; i < m_noNodes; i++)
		for(unsigned int j = 0; j < m_noNodes; j++)
		{
			double xx = m_newPheromoneMatrix->get(i,j);
			double yy = (*m_distanceMatrix)[i][j];
			(*m_heuristicMatrix)[i][j] = pow( xx , m_alpha) * pow(HEURISTIC(i,j),m_beta);
			(*m_heuristicMatrix)[j][i] = (*m_heuristicMatrix)[i][j];
		}	


}
void CAntColonySystem::updatePheromones()
{ 
	//evaporateAllPhero();
	globalPheromoneDeposit();

	calculateHeuristicMatrix();
}


void CAntColonySystem::localPheromoneUpdate(int ant, int step)
{
	int to = this->m_Ants[ant].getCity(step);
    int from = m_Ants[ant].getCity(step-1);
	double currentValue = this->m_newPheromoneMatrix->get(to, from);
	double updateValue =  (1.0 - xi) * currentValue + xi * tau0;

	this->m_newPheromoneMatrix->set(from, to , updateValue);
  	
	
	double dist = (*m_distanceMatrix)[to][from];
    (*m_heuristicMatrix)[from][to] = pow( updateValue,m_alpha)* pow(HEURISTIC(from,to),m_beta);
    (*m_heuristicMatrix)[to][from] = (*m_heuristicMatrix)[from][to];
}   
void CAntColonySystem::calculateNearestNeigbhor(unsigned int NUMOFANTS)
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
void CAntColonySystem::chooseClosestNext(std::vector<bool> &antsvisted , std::vector<int> &nntour)
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


void CAntColonySystem::initPheromones()
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
	
	
	//m_bestSoFarPath  = m_pLocalSearch->greedyPath(this->m_noNodes);

	//m_pLocalSearch->opt3(m_bestSoFarPath);

	//double distance = this->calculatePathLength(m_bestSoFarPath);
	//this->m_bestSoFarPathlength = (int)lenght;
	tau0 = 1.0 / ((double) m_noNodes * (double) distance);

	//tau0= 1.0 / ((double)lenght);
	for (unsigned int i=0; i<m_noNodes; i++)	
	{
		for (unsigned int j=0; j<m_noNodes; j++)
			this->m_newPheromoneMatrix->set(i, j , tau0);
	}

	for(unsigned int i=0;i<this->m_noNodes;i++)
		this->m_newPheromoneMatrix->set(i, i , 0.0);
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
				if (m_Ants[k].isCityVisited(i) == false && i!=c)	
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
			 if((maxHeuristic < choice && !(m_Ants[k].isCityVisited(j))) && c!=j)
			 {
				maxHeuristic = choice;
				maxHeuristicIdx = j;
			}
        }
		x = fRand(0, 1);	
        if(x < q0)  
		{
            m_Ants[k].setAntCity(step, maxHeuristicIdx);
            m_Ants[k].setCityVisited(maxHeuristicIdx);
        }
        else //if exploitation
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
	
		const std::vector<int>& tourvector = m_Ants[k].getAntsCityTour();
		m_Ants[k].setAntTourLength((int)calculatePathLength(tourvector));

		localPheromoneUpdate(k,m_noNodes);
	}
	updateBestSoFarPath();
}
