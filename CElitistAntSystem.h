#pragma once
#include "CAntSystem.h"

/*
A first improvement on the initial AS, called the elitist strategy for Ant System
(EAS), was introduced in Dorigo (1992) and Dorigo et al., (1991a, 1996). The idea is
to provide strong additional reinforcement to the arcs belonging to the best tour
found since the start of the algorithm; this tour is denoted as Tbs (best-so-far tour) in
the following. Note that this additional feedback to the best-so-far tour (which can
be viewed as additional pheromone deposited by an additional ant called best-so-far
ant) is another example of a daemon action of the ACO metaheuristic.
/*/

class CElitistAntSystem :
	public CAntSystem
{
	int mWeight; 
public:
	virtual void updatePheromones();
		
	CElitistAntSystem(int noNodes, int noAnts, int mWeight , MatrixArrayTypeInt *distanceMatrix);
	virtual ~CElitistAntSystem(void);
};

