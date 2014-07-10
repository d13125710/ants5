#include "StdAfx.h"
#include "MinMaxAntSystem.h"


MinMaxAntSystem::MinMaxAntSystem(int noNodes, int noAnts, int maxStagnationIterations , std::vector<std::vector<int> >  *Vdistance) 
	: AntSystem(noNodes , noAnts , *Vdistance)
{
	 m_stagnationIterations = 0;
     m_maxStagnationIterations = maxStagnationIterations;
     m_doOpt2 = false;
     m_doOpt3 = true;
	 m_opt = new LocalSearch(Vdistance);
	 m_bestdistance = std::numeric_limits<int>::max();
	 m_bestItterdistance = std::numeric_limits<int>::max();
	 restart_best_ant_tour_length =  std::numeric_limits<int>::max();
	
	 m_neg = 20;
	  if(m_neg > noNodes)
		 m_neg = noNodes-1;
	 lambda = 0.05;
	 reset_ant =false;
	 best_to_date_tour.resize(noNodes);  //best to date for even
	 restart_best_ant_path.resize(noNodes); ;
	 iteration_best_ant.resize(noNodes); ;

}
void MinMaxAntSystem::choose_closest_next(std::vector<bool> &antsvisted , std::vector<int> &nntour)
{
	int city, current_city, next_city, min_distance;
	next_city =this->m_noNodes;

	for(int phase = 1; phase < m_noNodes; phase++)
	{
		current_city = nntour[phase - 1];
		min_distance = std::numeric_limits<int>::max();; /* Search shortest edge */
		
			for (city = 0; city < m_noNodes; city++)
			{
				if (antsvisted[city])
				; /* city already visited */
				else
				{
				
					if ((*m_Vdistance)[current_city][city] < min_distance)
					{
						next_city = city;
						min_distance = (*m_Vdistance)[current_city][city];
					}
				}
			}
	
	
		nntour[phase] = next_city;
		antsvisted[next_city] = true;
	}

}

 double MinMaxAntSystem::node_branching(double l)
 {
	int i, m;
	double min, max, cutoff;
	double avg;

	std::vector<double> num_branches(m_noNodes);

	
	for (m = 0; m < m_neg; m++) 
	{
	   
	    min = (*m_Vpheromone)[m][nn_list[m][1]];
	    max = (*m_Vpheromone)[m][nn_list[m+1][1]];


	    for (i = 1; i < m_neg; i++) 
		{
			if ((*m_Vpheromone)[m][nn_list[m][i]] > max)
				max = (*m_Vpheromone)[m][nn_list[m][i]];
			if ((*m_Vpheromone)[m][nn_list[m][i]] < min)
				min = (*m_Vpheromone)[m][nn_list[m][i]];
	    }
	    cutoff = min + l * (max - min);

	    for (i = 1; i < m_neg; i++)
		{
			if ((*m_Vpheromone)[m][nn_list[m][i]] > cutoff)
				num_branches[m] += 1;
	    }
	}
	avg = 0.;
	for (m = 0; m < m_noNodes; m++) {
	    avg += num_branches[m];
	}
	
	return (avg / (double) (m_noNodes * 2));
}

void MinMaxAntSystem::calculateNearestNeigbhor(int NUMOFANTS)
{
	int i = 0;
	double *pHelpArray = new double[m_noNodes];
	double *pdistanceArray = new double[m_noNodes ];
	*(pdistanceArray) = std::numeric_limits<double>::max();

	nn_list.resize(m_noNodes);
		for (int i = 0; i < m_noNodes; ++i)
			nn_list[i].resize(NUMOFANTS);

	
	for (int node = 0; node < m_noNodes; node++) 
	{ 

			for (i = 0; i < m_noNodes; i++) 
			{
				*(pdistanceArray+i)=(*m_Vdistance)[node][i];
				*(pHelpArray+i)=  i;
			}
			double max = std::numeric_limits<double>::max() - 1; 
			*(pdistanceArray+node) =  std::numeric_limits<double>::max();  // set to a large value .. 
			this->m_opt->sort2(pdistanceArray, pHelpArray, 0, m_noNodes - 1);
			for (i = 0; i < NUMOFANTS; i++) 
			{
				nn_list[node][i] = *(pHelpArray+i);
			}
	}
	delete [] pHelpArray;
	delete [] pdistanceArray;
}

void MinMaxAntSystem::updateBestSoFarTour()
{
	
	    m_bestSoFarTourlength= (int) calculateTourLength(m_bestSoFarTour);

		
		if(this->m_iterations == 0) //first init
		{
			m_bestSoFarTour = getBestTourIrreration();
			m_bestSoFarTourlength= (int) calculateTourLength(m_bestSoFarTour);
			m_stagnationIterations =0;
		}

		int d = (int) getItterationBestTourLength() ;

		//best in itteration
		if (d < restart_best_ant_tour_length)
		{
				restart_best_ant_tour_length = d;
				restart_best_ant_path = getBestTourIrreration();
				m_stagnationIterations =0;
		}

		if(d < m_bestSoFarTourlength)
		{
            m_bestdistance = d;

			best_to_date_tour= getBestTourIrreration();  //best to date for even
			m_bestSoFarTour = getBestTourIrreration();
			m_stagnationIterations = 0;
        }
		else 
		{
            //copy the tour
			//iteration_best_ant
			m_bestItterdistance = d;
			m_iteration_best_ant = getBestTourIrreration();
			branching_factor = node_branching(lambda);
			double p_x = exp(log(0.05) / m_noNodes);
			trail_min = 1. * (1. - p_x) / (p_x * (double) ((this->m_noAnts + 1) / 2));
			trail_max = 1. / ((m_rho) * m_bestItterdistance);
			trail_0 = trail_max;
			trail_min = trail_max * trail_min;
			m_stagnationIterations ++;
	
        }
		
		

 }
void MinMaxAntSystem::initPheromones()
{
		std::vector<bool> visited(m_noNodes);
		std::vector<int> nntour(m_noNodes);
		calculateNearestNeigbhor(20);
		//calculate min max values inital
		int phase = 0;
		int rnd= (rand()%(visited.size()-1))+1;
		nntour[0] =rnd;
		visited[rnd] = true;
		choose_closest_next(visited, nntour);
		double distance = this->calculateTourLength( nntour);
		trail_max = 1. / ((m_rho) * distance);
		trail_min = trail_max / (2. * m_noNodes);
		init_pheromone_trails(trail_max);
}	
	
//There are variants in the selection of the ants allowed to update pheromones: the best to date ant, 
//or the best for current iteration, or the best after latest reset ant, or the best to date ant for even iterations, 
//and the best for iteration for odd iterations.

void MinMaxAntSystem::init_pheromone_trails(double initial_trail)
 {
	int i, j;

		for (i = 0; i < m_noNodes; i++)
		{
			for (j = 0; j <= i; j++) 
			{
				(*m_Vpheromone)[i][j] = initial_trail;
				(*m_Vpheromone)[j][i] = initial_trail;
			}
		}
}
void MinMaxAntSystem::global_update_pheromone(std::vector<int> &AntTour)
{
	int Q=1;

	double tour_length = calculateTourLength (AntTour); 
	for(int city = 1; city < AntTour.size(); city++)
	 {
	      int from = AntTour[city-1];
		  int to = AntTour[city];
 
		  // eq 14.2 / 14.3
		  (*m_Vpheromone)[from][to] += (Q/tour_length) ;
		  (*m_Vpheromone)[to][from] += (Q/tour_length) ; 
	
  }

}
void MinMaxAntSystem::check_pheromone_trail_limits()
{
	int i, j;

		for (i = 0; i < m_noNodes ; i++) 
		{
			for (j = 0; j < i; j++) 
			{
				if ( (*m_Vpheromone)[i][j] < trail_min)
				{
					 (*m_Vpheromone)[i][j] = trail_min;
					 (*m_Vpheromone)[j][i] = trail_min;
				} 
				else if ( (*m_Vpheromone)[i][j] > trail_max)
				{
					 (*m_Vpheromone)[i][j] = trail_max;
					 (*m_Vpheromone)[j][i] = trail_max;
				}
			}
		}
}


void MinMaxAntSystem::updatePheromones()
{
	if (this->m_iterations % 50 == 0) //even
	{
	  
		global_update_pheromone(best_to_date_tour);  //best to date for even
	} 
	else //odd
	{
	    if(reset_ant)
			global_update_pheromone(restart_best_ant_path);
	    else
			global_update_pheromone(getBestTourIrreration());
	}
	
    //end of ants walk
	Evaporation ();	
	check_pheromone_trail_limits();
			
	if ((m_iterations % 100) == 0 && m_stagnationIterations > 20) 
	{
		//create nn_list ..
		branching_factor = node_branching(lambda);
		double  branch_fac = 1.00001;
		if ( (branching_factor < branch_fac)	&& (m_iterations  > 150)) 
		{
			//cout << "INIT TRAILS!!!\n";
			restart_best_ant_tour_length =  std::numeric_limits<int>::max();
			init_pheromone_trails(trail_max);
			reset_ant = true;
			m_stagnationIterations =0;

		}

	}
	else
	{
		reset_ant = false;
	}

	this->calculateHeuristicMatrix();



    
}

void MinMaxAntSystem::localSearch()
{
        if(m_doOpt2){
            for(int i = 0; i < m_noAnts; i++){
			   m_opt->opt2(m_Ants[i].getCompleteTour());
            }
        }
        if(m_doOpt3){
            for(int i = 0; i < m_noAnts; i++){
                m_opt->opt3(m_Ants[i].getCompleteTour());
            }
        }

        updateBestSoFarTour();
        if(m_stagnationIterations > 0)
           m_stagnationIterations--;
    }



MinMaxAntSystem::~MinMaxAntSystem(void)
{
	if(m_opt)
		delete m_opt;
}
