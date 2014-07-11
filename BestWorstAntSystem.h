#pragma once
#include "cantsystem.h"
class CBestWorstAntSystem :
	public CAntSystem
{
	CAnt RestartBestAnt;
	int worstIndex;
	int Restart;

public:
	CBestWorstAntSystem(int noNodes, int noAnts , MatrixArrayTypeInt *distanceMatrix);
	virtual ~CBestWorstAntSystem(void);

	
public:
	virtual void updatePheromones();
	void initPheromoneTrails(double initialValue) const;
	void bwas_worst_ant_update( const CAnt &bestAnt, const CAnt &worstAnt);
	long distance_between_ants( const CAnt &bestAnt, const CAnt &worstAnt2);
	void bwas_pheromone_mutation(  );

};

