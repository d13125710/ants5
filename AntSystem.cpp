#include "StdAfx.h"
#include "AntSystem.h"


AntSystem::AntSystem(int noNodes, int noAnts)
{
	Create(noNodes, noAnts);
}

AntSystem::AntSystem(int noNodes, int noAnts, MatrixArrayTypeInt &distanceMatrix)
{
	m_Vdistance = &distanceMatrix;
	Create(noNodes, noAnts);
	
	      
}
void AntSystem::Create(int noNodes, int noAnts)
{
	if(noAnts > noNodes)
		noAnts = noNodes-1;
	
	m_opt = new LocalSearch(m_Vdistance);


	this->m_noNodes = noNodes;
    this->m_noAnts = noAnts;
    this->m_alpha = 1.0;
    this->m_beta = 5.0;
    this->m_rho = 0.5;
    this->m_iterations = 0;

	//m_doOpt3
	m_doOpt2 = false;
	m_doOpt3 = false;
	m_doGreed = true;

	m_Vpheromone = new std::vector<std::vector<double> >(noNodes , noNodes);
	m_Vchoice= new std::vector<std::vector<double> >(noNodes , noNodes);
	m_bestSoFarTour.resize(noNodes);
	m_Ants.empty();
	m_Ants.resize(noAnts);
	for(int i = 0; i < noAnts; i++)
	{
		Ant pant(noNodes);
		m_Ants[i] = pant;
	}

	//start ramdom
	//srand((unsigned int)time(0));
	//std::random_device random_dev;
	
	//std::mt19937 mt(1729);
	
	

	m_bestIrrTour.resize(m_noNodes);
	m_Randomnseq.resize(m_noNodes);
			
	for (int i=0; i<m_noNodes; i++)
		m_Randomnseq[i]=i;
	std::random_shuffle( m_Randomnseq.begin() , m_Randomnseq.end() );
	//m_generator((random_dev());
	t_prob.resize(noNodes);
	strength.resize(noNodes+1);

	m_bestSoFarTourlength = std::numeric_limits<int>::max();

	
	 
      
}

void  AntSystem::updatePheromones()
{
        Evaporation();
        for(int k = 0; k < m_noAnts; k++)
            adjust_pheremone(k);
		//dont think i need this on irreration
        calculateHeuristicMatrix();
}

void AntSystem::initMatrix(std::vector<std::vector<int> > &distanceMatrix)
{
    int i,j;
	this->m_Vdistance = &distanceMatrix;
     
	for(i = 0; i < m_noNodes; i++)
            for(j = 0; j < m_noNodes; j++)
			{
                (*m_Vpheromone)[i][j] = 0.0;
                (*m_Vchoice)[i][j] = 0.0;
            }
}
 void AntSystem::initPheromones()
{
	for (int i=0; i<m_noNodes; i++)	{
		for (int j=0; j<m_noNodes; j++)
			(*m_Vpheromone)[i][j] =  0.1 * rand() / (double)RAND_MAX;
	}
       
}
 void AntSystem::initAnts()
 {
     for(unsigned int i = 0; i < m_Ants.size(); i++)
		m_Ants[i].clearAntsMemory();

	
	shuffle(m_Randomnseq.begin(), m_Randomnseq.end(), m_generator);

 }

void AntSystem::adjust_pheremone(int antPos)
{
	Ant &pant=m_Ants[antPos];

	double tour_length = 1 / calculateTourLength(pant.getCompleteTour()); 
	for(int city = 1; city < m_noNodes; city++)
	 {
	      
		int from = m_Ants[antPos].getTour(city);
		int to   = m_Ants[antPos].getTour(city+1);
		// eq 14.2 / 14.3
		(*m_Vpheromone)[from][to] += tour_length ;
		(*m_Vpheromone)[to][from] += tour_length ; 
	
  }

}

void AntSystem::decisionRule(int k, int step)
{
        int c = m_Ants[k].getTour(step-1); 
 
		for (int i=0; i<m_noNodes; i++) 
		{
			t_prob[i] = 0.0;
			strength[i] = 0;	
			if (m_Ants[k].isVisited(i) == false)	
			{
				t_prob[i]= (*m_Vchoice)[c][i];
						
			}
		}

		for (int z =0; z < m_noNodes; z++)
			strength[z+1] = t_prob[z] + strength[z];

		
		
		
		long _rnd = rand();

		if (_rnd == 0)
			_rnd++;

		double x = strength[m_noNodes] * (double)_rnd / RAND_MAX;	
				
		int y = 0;
		while (!((strength[y] <= x) && (x <= strength[y+1])))
			y++;
			
		m_Ants[k].setTour(step, y);
		m_Ants[k].setVisited(y);

}


void AntSystem::constructSolutions()
{
	initAnts();
	   
	for (int k = 0; k < m_noAnts; k++ )
	{
		m_Ants[k].setTour(0,m_Randomnseq[k]);
		m_Ants[k].setVisited(m_Randomnseq[k]);

	}

	for(int step = 1 ; step < m_noNodes; step++)
	{
		for(int k = 0; k < m_noAnts; k++)
			decisionRule(k,step);
	}
     
	for(int k = 0; k < m_noAnts; k++){
			
		int tourstart=m_Ants[k].getTour(0);
		m_Ants[k].setTour(m_noNodes,tourstart);
		const std::vector<int>& tourvector = m_Ants[k].getCompleteTour();
		m_Ants[k].settourLength((int)calculateTourLength(tourvector));
	}
	updateBestSoFarTour();
}

std::vector<int> &AntSystem::getBestTourIrreration()
{
    
    int bestTourLength = std::numeric_limits<int>::max();

	for(int i = 0; i < m_Ants.size() ; i++)
	{
        if(m_Ants[i].getAntTourLength() < bestTourLength)
		{
             bestTourLength = m_Ants[i].getAntTourLength();
             for(int j = 0; j < m_noNodes; j++)
			 {
                    m_bestIrrTour[j] = m_Ants[i].getTour(j);
			 }
        }
    }
    return m_bestIrrTour;
}
void AntSystem::Evaporation()
{
	int from, to;
	for (from=0; from < m_noNodes; from++) {
		for (to=0; to<from; to++) {
			(*m_Vpheromone)[from][to] *= (1-m_rho);
			(*m_Vpheromone)[to][from] *= (1-m_rho);

		}
	}
}

double AntSystem::calculateTourLength(const std::vector<int> &currentTour)
{
    int  start, end;
	double d = 0;
	for ( int i=1; i<m_noNodes; i++)
	{
		start = currentTour[i-1];
		end = currentTour[i];
		d += (*m_Vdistance)[start][end];	
	}

	d += (*m_Vdistance)[currentTour[m_noNodes-1]][currentTour[0]];
	return d;
}
void AntSystem::calculateHeuristicMatrix()
{
       for(int i = 0; i < m_noNodes; i++)
          for(int j = 0; j < m_noNodes; j++)
		  {
                double xx = (*m_Vpheromone)[i][j];
				double yy = (*m_Vdistance)[i][j];

				(*m_Vchoice)[i][j] = pow( xx , m_alpha) / pow(yy,m_beta);
		  }	
}
void AntSystem::updateBestSoFarTour()
{
		int distance= m_bestSoFarTourlength ;;//(int) calculateTourLength(m_bestSoFarTour);
		int bestitter = getItterationBestTourLength() ;
		if(distance == 0)
		{
			m_bestSoFarTour = getBestTourIrreration();
			m_bestSoFarTourlength = getItterationBestTourLength() ;

		}
        else if(bestitter < distance)
		{
            m_bestSoFarTourlength = bestitter;
			m_bestSoFarTour = getBestTourIrreration();
        }
}


AntSystem::~AntSystem(void)
{
	delete m_Vpheromone;
	delete m_Vchoice;
	delete m_opt;
}
