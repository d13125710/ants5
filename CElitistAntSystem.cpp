#include "StdAfx.h"
#include "CElitistAntSystem.h"


//get the best distance and adjust phero with that distance
//reported, that the best value for mWeight is between four and six.
CElitistAntSystem::CElitistAntSystem(int noNodes, int noAnts, int weight , MatrixArrayTypeInt *distanceMatrix) 
	: CAntSystem(noNodes , noAnts , distanceMatrix)
{
	this->mWeight = weight; 
}
//The pheromone update is a little bit different: on each iteration, the best to date ant deposits an additional quantity of pheromone on paths it traveled:
void CElitistAntSystem::updatePheromones()
{
	evaporateAllPhero();
	
	//add best tour stuff
	double d_tau = (double)mWeight/ m_bestSoFarPathlength; 
	for(int city = 1; city <m_bestSoFarPath.size(); city++)
	{
		int from = m_bestSoFarPath[city-1];
		int to = m_bestSoFarPath[city];
 		// eq 14.2 / 14.3
		this->m_newPheromoneMatrix->add(from , to , d_tau);

		(*m_pheromoneMatrix)[from][to]+= d_tau;  //sermertic array
	    (*m_pheromoneMatrix)[to][from]=  (*m_pheromoneMatrix)[from][to]; 
	
	}
	//for(std::vector<CAnt>::iterator it=m_Ants.begin(); it!=m_Ants.end();it++) {
	//CAnt &ant = (*it);
////	ant.
//}
	

	for(unsigned int k = 0; k < m_noAnts; k++)
		adjustPheremone(k);

	calculateHeuristicMatrix();
}

CElitistAntSystem::~CElitistAntSystem(void)
{
}
