#include "StdAfx.h"
#include "Chromo1.h"


Chromo1::Chromo1(int noGenes, MatrixArrayTypeInt *distance)
{
	this->m_noGenes= noGenes;
	m_dist = distance;

	tour.resize(noGenes);

	fitness = 0;

	for(int i=0;i<tour.size();i++)
		tour[i]=-1;

	//start ramdom
	
	//m_generator((random_dev());
}


Chromo1::~Chromo1(void)
{
}
