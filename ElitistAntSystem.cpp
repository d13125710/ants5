#include "StdAfx.h"
#include "ElitistAntSystem.h"


//get the best distance and adjust phero with that distance
//reported, that the best value for mWeight is between four and six.
ElitistAntSystem::ElitistAntSystem(int noNodes, int noAnts, int weight , MatrixArrayTypeInt *distanceMatrix) 
	: AntSystem(noNodes , noAnts , *distanceMatrix)
{
	this->mWeight = weight;
}
//The pheromone update is a little bit different: on each iteration, the best to date ant deposits an additional quantity of pheromone on paths it traveled:
void ElitistAntSystem::updatePheromones()
{
	
	
	//add best tour stuff
	double d_tau = (double)mWeight/ m_bestSoFarTourlength; 
	for(int city = 1; city <m_bestSoFarTour.size(); city++)
	{
		int from = m_bestSoFarTour[city-1];
		int to = m_bestSoFarTour[city];
 		// eq 14.2 / 14.3
		(*m_Vpheromone)[from][to]+= d_tau;  //sermertic array
	    (*m_Vpheromone)[to][from]=  (*m_Vpheromone)[from][to]; 
	
	}

	//call the base class update first
	AntSystem::updatePheromones();
}

ElitistAntSystem::~ElitistAntSystem(void)
{
}
