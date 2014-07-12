#include "StdAfx.h"
#include "BestWorstAntSystem.h"



CBestWorstAntSystem::CBestWorstAntSystem(int noNodes, int noAnts,MatrixArrayTypeInt *distanceMatrix) 
	: CAntSystem(noNodes , noAnts , distanceMatrix)
{
	Restart =0;
}
/*    
      FUNCTION:      uses additional evaporation on the arcs of iteration worst
                     ant that are not shared with the global best ant
      INPUT:         pointer to the worst (a1) and the best (a2) ant
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones on some arcs undergo additional evaporation
*/
void CBestWorstAntSystem::bwas_worst_ant_update( const CAnt &worseA, const CAnt &bestA)
{  
    unsigned int   i, j, h, pos, pred;
    unsigned int   distance;
    
	std::vector<int> pos2(m_noNodes);
	
	

    for ( i = 0 ; i < m_noNodes ; i++ ) 
	{
		pos2[bestA.getCity(i)] = i;
    }
 
    distance = 0;
    for ( i = 0 ; i < m_noNodes -1  ; i++ ) 
	{
		j = worseA.getCity(i);
		h = worseA.getCity(i+1);
		pos = pos2[j];
		if (pos - 1 < 0)
			pred = m_noNodes - 1;
		else 
			pred = pos - 1;
		if (bestA.getCity(pos+1) == h)
			; /* do nothing, edge is common with a2 (best solution found so far) */
		else if (bestA.getCity(pred) == h)
			; /* do nothing, edge is common with a2 (best solution found so far) */
		else {   /* edge (j,h) does not occur in ant a2 */       
			m_newPheromoneMatrix->evaporate(j,h);
			
			//(*m_pheromoneMatrix)[j][h] = (1 - m_rho) * (*m_pheromoneMatrix)[j][h];
			//(*m_pheromoneMatrix)[h][j] = (1 - m_rho) * (*m_pheromoneMatrix)[h][j];
		}
    }

}
long CBestWorstAntSystem::distance_between_ants( const CAnt &bestAnt, const CAnt &worstAnt)
{
	unsigned int j, h, pos, pred;
	std::vector<int> pos2(m_noNodes);
	

	for (unsigned int i = 0 ; i < m_noNodes ; i++ ) 
	{
		pos2[worstAnt.getCity(i)] = i;
	}

	long distance = 0;
	for (unsigned int i = 0; i < m_noNodes -1; i++) {
		j = bestAnt.getCity(i);
		h = bestAnt.getCity(i + 1);
		pos = pos2[j];
		if (pos - 1 < 0)
			pred = m_noNodes - 1;
		else
			pred = pos - 1;
		if (worstAnt.getCity(i + 1)== h)
			; /* do nothing, edge is common with best solution found so far */
		else if (worstAnt.getCity(pred) == h)
			; /* do nothing, edge is common with best solution found so far */
		else { /* edge (j,h) does not occur in ant a2 */
			distance++;
		}
	}
	return distance;
}
void CBestWorstAntSystem::initPheromoneTrails(double initialValue) const
{
	unsigned int i, j;

	for (i = 0; i < m_noNodes; i++)
	{
		for (j = 0; j <= i; j++) 
		{
			(*m_pheromoneMatrix)[i][j] = initialValue;
			(*m_pheromoneMatrix)[j][i] = initialValue;
		}
	}
}

//The pheromone update is a little bit different: on each iteration, the best to date ant deposits an additional quantity of pheromone on paths it traveled:
void CBestWorstAntSystem::updatePheromones()
{
	//get the best at so far

	//this->evaporateAllPhero();
	

	//global_update_pheromone( best_so_far_ant );
	double d_tau = (double)1/ m_bestSoFarPathlength; 
	for(int city = 1; city <m_bestSoFarPath.size(); city++)
	{
		int from = m_bestSoFarPath[city-1];
		int to = m_bestSoFarPath[city];
		// eq 14.2 / 14.3
		(*m_pheromoneMatrix)[from][to]+= d_tau;  //sermertic array
		(*m_pheromoneMatrix)[to][from]=  (*m_pheromoneMatrix)[from][to]; 

	}

	//iteration_worst_ant = find_worst();
	int worstDistance = m_bestSoFarPathlength;
	
	for(unsigned int i =0; i< this->m_noNodes ; i++)
	{
		if(this->m_Ants[i].getAntTourLength() > worstDistance)
		{
			 worstIndex= i;
			 worstDistance = this->m_Ants[i].getAntTourLength();

		}
	}
	//bwas_worst_ant_update
	bwas_worst_ant_update( m_Ants[worstIndex], GetBestAnt() );
	
	
	//mutate every 150 itter
	//if(irreationSinceLastBest > m_noNodes)
	//{
		long distance_best_worst = distance_between_ants(  GetBestAnt() ,m_Ants[worstIndex]);
		printf("distance_best_worst %ld, tour length worst %ld\n",distance_best_worst, m_Ants[worstIndex].getAntTourLength()); 
		if ( distance_best_worst < (long int) (0.05 * (double) m_noNodes) ) 
		{
			m_bestSoFarPathlength =std::numeric_limits<int>::max();;
			RestartBestAnt.setAntTourLength(std::numeric_limits<int>::max());
			initPheromones(); //random
			printf("init pheromone trails with iteration %ld\n" , this->m_iterations);
		}
		else 
			bwas_pheromone_mutation();
		irreationSinceLastBest=0;
	//}
	
	calculateHeuristicMatrix();
}
void CBestWorstAntSystem::bwas_pheromone_mutation( void )
/*    
      FUNCTION: implements the pheromone mutation in Best-Worst Ant System
      INPUT:    none  
      OUTPUT:   none
*/
{
   unsigned int     i, j, k;
    unsigned int     num_mutations;
    double       avg_trail = 0.0, mutation_strength = 0.0, mutation_rate = 0.3;

     /* compute average pheromone trail on edges of global best solution */
    for ( i = 0 ; i < m_noNodes-1 ; i++ ) 
	{
		avg_trail +=  (*m_pheromoneMatrix)[m_bestSoFarPath[i]][m_bestSoFarPath[i+1]];
    }
    avg_trail /= (double) m_noNodes;
  
    /* determine mutation strength of pheromone matrix */ 
    /* FIXME: we add a small value to the denominator to avoid any
       potential division by zero. This may not be fully correct
       according to the original BWAS. */
    //if ( max_time > 0.1 )
	 mutation_strength =  avg_trail / fRand(0.1,0.9);
   // else if ( max_tours > 100 )
//	mutation_strength = 4. * avg_trail * (iteration - restart_iteration) 
 //           / (max_tours - restart_iteration + 1);
  //  else
//	printf("apparently no termination condition applied!!\n");

    /* finally use fast version of matrix mutation */
    mutation_rate =  0.1 * rand() / (double)RAND_MAX;;
    num_mutations = static_cast<long>(m_noNodes * mutation_rate / 2);   
    /* / 2 because of adjustment for symmetry of pheromone trails */
 
   // if ( restart_iteration < 2 )
	//num_mutations = 0; 

    for ( i = 0 ; i < num_mutations ; i++ ) 
	{
		j =   (long int) (fRand(0,1) * (double) m_noNodes);
		k =   (long int) (fRand(0,1) * (double) m_noNodes);
		if ( fRand(0,1) < 0.5 ) {
			(*m_pheromoneMatrix)[j][k] += mutation_strength;
			(*m_pheromoneMatrix)[k][j] = (*m_pheromoneMatrix)[j][k];
		}
		else {
			(*m_pheromoneMatrix)[j][k] -= mutation_strength;
			if ( (*m_pheromoneMatrix)[j][k] <= 0.0 ) {
				(*m_pheromoneMatrix)[j][k] = 0.00000000000000000000000000000001;
			}
			(*m_pheromoneMatrix)[k][j] = (*m_pheromoneMatrix)[j][k]; 
		}
    }
}




CBestWorstAntSystem::~CBestWorstAntSystem(void)
{
}
