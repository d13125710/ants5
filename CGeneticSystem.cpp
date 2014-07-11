#include "StdAfx.h"
#include "CGeneticSystem.h"


CGeneticSystem::CGeneticSystem(void)
{

}


CGeneticSystem::~CGeneticSystem(void)
{
	delete m_pLocalSearch;
}

CGeneticSystem::CGeneticSystem(int populationSize, int noNodes, double mutationRatio, bool doOpt2, bool doOpt3 , MatrixArrayTypeInt *d)
{
	this->m_populationSize = populationSize;
	this->m_noCitys = noNodes;
	this->m_mutationRatio = mutationRatio;
	this->m_bOpt2Search = doOpt2;
	this->m_bOpt3Search = doOpt3;
	m_ChromoPopulation.resize(populationSize);
	firsttime =1;
	m_distMatrix = d;
	for(int i = 0; i < populationSize; i++)
	{
		CChromo p(noNodes , m_distMatrix);
		m_ChromoPopulation[i] = p;
	}
	m_fitness.resize(populationSize); 
	m_bestSoFarPath.resize(noNodes);
	for(int i = 0; i < noNodes; i++)
		m_bestSoFarPath[i] = i;
	m_currentIteration = 0;
	m_pLocalSearch = new CLocalSearch(m_distMatrix);
	}
    
void CGeneticSystem::setIteration(int iter)
{
	m_currentIteration = iter;
}
    
int CGeneticSystem::getIteration() const
{
	return m_currentIteration;
}
    
void CGeneticSystem::initPopulation()
{
	for(int i = 0; i < m_populationSize; i++)
	{
		m_ChromoPopulation[i].shuffle();
	}

}

void CGeneticSystem::updateBestSoFarPath()
{
   
	int dist = computePathLength(m_bestSoFarPath);

	int irrdist = getBestPathLength();
	if(dist > irrdist || firsttime == 1)
	{
		  m_bestSoFarPath =  m_ChromoPopulation[0].getPath();
	}


//	for(int i = 0; i < m_populationSize; i++)
//	{
//		if(computePathLength(m_bestSoFarPath) > m_ChromoPopulation[i].getDistance() || this->firsttime ==1)
//		{
//			firsttime = 0;
//            for(int j = 0; j < m_noCitys; j++)
//                m_bestSoFarPath[j] = m_ChromoPopulation[i].getGene(j);
  //      }

//    }
}
std::vector<int> & CGeneticSystem::getBestPath()
{
	SortPopulation(m_ChromoPopulation , true);
	return m_ChromoPopulation[0].getPath();
}
int CGeneticSystem::getBestPathLength()
{
	SortPopulation(m_ChromoPopulation , true);
 	return m_ChromoPopulation[0].getDistance();
}
int CGeneticSystem::computePathLength(const std::vector<int> &Path)
{
	int len = 0;
	for(int i = 0; i < m_noCitys-1; i++){
		len+=(*m_distMatrix)[Path[i]][Path[i+1]];
	}
	len+=(*m_distMatrix)[Path[m_noCitys-1]][Path[0]];
	return len;
}
// Selects candidate tour for crossover
int CGeneticSystem::tournamentSelection()
{
	double sumFitness = 0.0;
	int idx1 = 0;
    for(int i = 0; i < m_populationSize; i++)
        sumFitness += m_fitness[i];
	double _rnd = rand();
	if(_rnd==0)
		_rnd++;
	double prob = ((double)_rnd / RAND_MAX) * sumFitness;
		
    double sum = m_fitness[0];
    while(sum < prob)
	{
        idx1++;
        sum+=m_fitness[idx1];
    }
        
	return idx1;	
        
}
void CGeneticSystem::stepGeneration2()
{
	std::vector<CChromo> newPopulation(m_populationSize*2);
	for(int i = 0; i < m_populationSize*2; i++)
	{
		CChromo p(m_noCitys , m_distMatrix);
		newPopulation[i] = p;
	}

	int newPopulationSize = 0;
	SortPopulation(m_ChromoPopulation , false);
	computeFitness();

	for(int i = 0; i < m_populationSize; i++)
	{
		for(int j = 0; j < m_noCitys; j++)
		{
			int test = m_ChromoPopulation[i].getGene(j);
			newPopulation[i].setGene(j,test );
		}
		newPopulationSize++;
	}
     

	while (newPopulationSize < 2*m_populationSize) 
	{
		int idx1 = tournamentSelection();
		int idx2 = tournamentSelection();
		CChromo offspring = m_ChromoPopulation[idx1].CrossOver2(&m_ChromoPopulation[idx2]);
		newPopulation[newPopulationSize] = offspring;
		newPopulationSize++;
	}
  	mutatePopulation(newPopulation);
	SortPopulation(newPopulation , true);
	for(int i = 0; i < m_populationSize; i++)
	{
		m_ChromoPopulation[i] = newPopulation[i];
	}
	SortPopulation(m_ChromoPopulation,true);
	updateBestSoFarPath();
}



void CGeneticSystem::stepGeneration(){
	//stepGeneration2();
	//return;
	std::vector<CChromo> newPopulation(m_populationSize*2);
	for(int i = 0; i < m_populationSize*2; i++)
	{
		CChromo p(m_noCitys , m_distMatrix);
		newPopulation[i] = p;
	}
    int newPopulationSize = 0;
    SortPopulation(m_ChromoPopulation,false);
    computeFitness();
    for(int i = 0; i < m_populationSize; i++)
	{
	  for(int j = 0; j < m_noCitys; j++)
		{
			int test = m_ChromoPopulation[i].getGene(j);
            newPopulation[i].setGene(j,test );
        }
        newPopulationSize++;
    }
    while(newPopulationSize < 2*m_populationSize)
	{
        int idx1 =0; //tournamentSelection();
		int idx2 =0; //tournamentSelection();
		while(idx1 == idx2)
		{
			idx2= tournamentSelection();
			idx1= tournamentSelection();
		}
		CChromo &pfather=m_ChromoPopulation[idx2];
		CChromo &pMother=m_ChromoPopulation[idx1];
        CChromo p_offspring1(m_noCitys,m_distMatrix) , p_offspring2(m_noCitys,m_distMatrix);
		pMother.crossover(&pfather, &p_offspring1, &p_offspring2);
		newPopulation[newPopulationSize] = p_offspring1;
        newPopulationSize++;
		if(newPopulationSize >= newPopulation.size())
            break;
        newPopulation[newPopulationSize] = p_offspring2;
        newPopulationSize++;

    }
    mutatePopulation(newPopulation);
    SortPopulation(newPopulation , true);  //ass
    for(int i = 0; i < m_populationSize-2; i++) //keep last best
	{
		m_ChromoPopulation[i] = newPopulation[i];
	}
    SortPopulation(m_ChromoPopulation , true);
    updateBestSoFarPath();
}
void CGeneticSystem::computeFitness()
{
    double selectionPressure = 2;
    for(int i = 0; i < m_populationSize; i++)
	{
        double fit = (double) 1*(i+1)/ (double) (m_populationSize-1);
        m_fitness[i] = fit;
    }
}
	

void CGeneticSystem::mutatePopulation(std::vector<CChromo> &popN)
{
    for(int i = 0; i < popN.size(); i++)
	{
		double rnd =((double) rand() / (RAND_MAX)) ;
		if(rnd < m_mutationRatio)
			popN[i].mutate();
    }
}

void CGeneticSystem::SortPopulation(std::vector<CChromo> &pop , bool Direction)
{
	bool done = false;
	while(!done){
		done = true;
		for(int i = 0; i < pop.size()-1; i++)
		{
			int from = (Direction) ? pop[i].getDistance() : pop[i+1].getDistance();
			int to = (Direction) ? pop[i+1].getDistance() : pop[i].getDistance();
			
			if(from > to)  //true
			{
				CChromo swap =pop[i];
				//swap.
				pop[i] = pop[i+1];
				pop[i+1] = swap;
				done = false;
			}
		}
	}
	
}
const std::vector<int> & CGeneticSystem::getBestSoFarPath()
{
	return m_bestSoFarPath;
}
int CGeneticSystem::getBestSoFarPathLength()
{
    return computePathLength(m_bestSoFarPath);
}
void CGeneticSystem::localSearch(){
    if(m_bOpt2Search)
	{
        for(int i = 0; i < m_populationSize; i++)
		{
            m_pLocalSearch->opt2(m_ChromoPopulation[i].getPath());
        }
    }
    if(m_bOpt3Search)
	{
        for(int i = 0; i < m_populationSize; i++)
		{
            m_pLocalSearch->opt3(m_ChromoPopulation[i].getPath());
        }
    }
    updateBestSoFarPath();
}

