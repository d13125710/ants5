#include "StdAfx.h"
#include "CAntSystem.h"


CAntSystem::CAntSystem(int noNodes, int noAnts)
{
	Create(noNodes, noAnts);
}

CAntSystem::CAntSystem(int noNodes, int noAnts, MatrixArrayTypeInt *distanceMatrix)
{
	m_distanceMatrix = distanceMatrix;
	Create(noNodes, noAnts);
	
	      
}
void CAntSystem::Create(int noNodes, int noAnts)
{
	if(noAnts > noNodes)
		noAnts = noNodes-1;
	irreationSinceLastBest=0;
	m_pLocalSearch = new CLocalSearch(m_distanceMatrix);


	this->m_noNodes = noNodes;
    this->m_noAnts = noAnts;
    this->m_alpha = 1.0;
    this->m_beta = 2.0;
    this->m_rho = 0.5;
    this->m_iterations = 0;

	//m_doOpt3
	m_LocalSearchOpt2 = false;
	m_LocalSearchOpt3 = false;
	m_LocalSearchGreed = true;

	m_pheromoneMatrix = new std::vector<std::vector<double> >(noNodes , noNodes);
	m_heuristicMatrix= new std::vector<std::vector<double> >(noNodes , noNodes);
	m_bestSoFarPath.resize(noNodes);
	m_Ants.empty();
	m_Ants.resize(noAnts);
	for(int i = 0; i < noAnts; i++)
	{
		CAnt newAnt(noNodes);
		m_Ants[i] = newAnt;
	}
	m_bestIrrPath.resize(m_noNodes);
	m_randomPath.resize(m_noNodes);
			
	for (int i=0; i<m_noNodes; i++)
		m_randomPath[i]=i;
	std::random_shuffle( m_randomPath.begin() , m_randomPath.end() );
	//m_generator((random_dev());
	t_prob.resize(noNodes);
	m_strength.resize(noNodes+1);

	m_bestSoFarPathlength = std::numeric_limits<int>::max();

	
	 
      
}

void  CAntSystem::updatePheromones()
{
    evaporateAllPhero();
    for(int k = 0; k < m_noAnts; k++)
        adjustPheremone(k);
	//dont think i need this on irreration
    calculateHeuristicMatrix();
}

void CAntSystem::initMatrix(std::vector<std::vector<int> > &distanceMatrix)
{
    int i,j;
	this->m_distanceMatrix = &distanceMatrix;
     
	for(i = 0; i < m_noNodes; i++)
            for(j = 0; j < m_noNodes; j++)
			{
                (*m_pheromoneMatrix)[i][j] = 0.0;
                (*m_heuristicMatrix)[i][j] = 0.0;
            }
}
 void CAntSystem::initPheromones()
{
	for (int i=0; i<m_noNodes; i++)	{
		for (int j=0; j<m_noNodes; j++)
			(*m_pheromoneMatrix)[i][j] =  0.1 * rand() / (double)RAND_MAX;
	}
       
}
 void CAntSystem::initAnts()
 {
     for(unsigned int i = 0; i < m_Ants.size(); i++)
		m_Ants[i].ReinitiseAnt();

	
//	shuffle(m_Randomnseq.begin(), m_Randomnseq.end(), m_generator);

 }

void CAntSystem::adjustPheremone(int antPos)
{
	CAnt &pant=m_Ants[antPos];

	double tour_length = 1 / calculatePathLength(pant.getAntsCityTour()); 

	
	for(int city = 1; city < m_noNodes; city++)
	 {
	      
		int from = m_Ants[antPos].getCity(city);
		int to   = m_Ants[antPos].getCity(city+1);
		// eq 14.2 / 14.3
		(*m_pheromoneMatrix)[from][to] += tour_length ;
		(*m_pheromoneMatrix)[to][from] += tour_length ; 
	
  }

}

void CAntSystem::decisionRule(int k, int step)
{
        int c = m_Ants[k].getCity(step-1); 
 
		for (int i=0; i<m_noNodes; i++) 
		{
			t_prob[i] = 0.0;
			m_strength[i] = 0;	
			if (m_Ants[k].isCityVisited(i) == false)	
			{
				t_prob[i]= (*m_heuristicMatrix)[c][i];
						
			}
		}

		for (int z =0; z < m_noNodes; z++)
			m_strength[z+1] = t_prob[z] + m_strength[z];

		
		double x = fRand(0,  m_strength[m_noNodes]);
		
		
		
				
		int y = 0;
		while (!((m_strength[y] <= x) && (x <= m_strength[y+1])))
			y++;
			
		m_Ants[k].setAntCity(step, y);
		m_Ants[k].setCityVisited(y);

}


void CAntSystem::constructSolutions()
{
	initAnts();
	   
	for (int k = 0; k < m_noAnts; k++ )
	{
		m_Ants[k].setAntCity(0,m_randomPath[k]);
		m_Ants[k].setCityVisited(m_randomPath[k]);

	}

	for(int step = 1 ; step < m_noNodes; step++)
	{
		for(int k = 0; k < m_noAnts; k++)
			decisionRule(k,step);
	}
     
	for(int k = 0; k < m_noAnts; k++){
			
		int tourstart=m_Ants[k].getCity(0);
		m_Ants[k].setAntCity(m_noNodes,tourstart);
		const std::vector<int>& tourvector = m_Ants[k].getAntsCityTour();
		m_Ants[k].setAntTourLength((int)calculatePathLength(tourvector));
	}
	updateBestSoFarPath();
}

std::vector<int> &CAntSystem::getBestTourIrreration()
{
    
    int bestTourLength = std::numeric_limits<int>::max();

	for(int i = 0; i < m_Ants.size() ; i++)
	{
        if(m_Ants[i].getAntTourLength() < bestTourLength)
		{
             bestTourLength = m_Ants[i].getAntTourLength();
             BestAntIndex = i;
			 for(int j = 0; j < m_noNodes; j++)
			 {
                    m_bestIrrPath[j] = m_Ants[i].getCity(j);
					
			 }
        }
    }
    return m_bestIrrPath;
}
void CAntSystem::evaporateAllPhero()
{
	int from, to;
	for (from=0; from < m_noNodes; from++) {
		for (to=0; to<=from; to++) {
			(*m_pheromoneMatrix)[from][to] *= (1-m_rho);
			(*m_pheromoneMatrix)[to][from] *= (1-m_rho);

		}
	}
}

double CAntSystem::calculatePathLength(const std::vector<int> &currentPath) const
{
    int  start, end;
	double d = 0;
	for ( int i=1; i<m_noNodes; i++)
	{
		start = currentPath[i-1];
		end = currentPath[i];
		d += (*m_distanceMatrix)[start][end];	
	}

	d += (*m_distanceMatrix)[currentPath[m_noNodes-1]][currentPath[0]];
	return d;
}
void CAntSystem::calculateHeuristicMatrix()
{
       for(int i = 0; i < m_noNodes; i++)
          for(int j = 0; j < m_noNodes; j++)
		  {
                double xx = (*m_pheromoneMatrix)[i][j];
				double yy = (*m_distanceMatrix)[i][j];

				(*m_heuristicMatrix)[i][j] = pow( xx , m_alpha) / pow(yy,m_beta);
		  }	
}
void CAntSystem::updateBestSoFarPath()
{
		int distance= m_bestSoFarPathlength ;;//(int) calculateTourLength(m_bestSoFarTour);
		int bestitter = getItterationBestPathLength() ;
		if(distance == 0)
		{
			m_bestSoFarPath = getBestTourIrreration();
			m_bestSoFarPathlength = getItterationBestPathLength() ;
			irreationSinceLastBest=0;
		}
        else if(bestitter < distance)
		{
			irreationSinceLastBest=0;
			m_bestSoFarPathlength = bestitter;
			m_bestSoFarPath = getBestTourIrreration();
        }
		else
		{
			irreationSinceLastBest++;
		}
}


CAntSystem::~CAntSystem(void)
{
	delete m_pheromoneMatrix;
	delete m_heuristicMatrix;
	delete m_pLocalSearch;
}
