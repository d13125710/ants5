#include "StdAfx.h"
#include "GeneticSystem.h"


GeneticSystem::GeneticSystem(void)
{

}


GeneticSystem::~GeneticSystem(void)
{
	delete opt;
}

 GeneticSystem::GeneticSystem(int populationSize, int noNodes, double mutationRatio, bool doOpt2, bool doOpt3 , MatrixArrayTypeInt *d)
 {
        this->populationSize = populationSize;
        this->noNodes = noNodes;
        this->mutationRatio = mutationRatio;
        this->doOpt2 = doOpt2;
        this->doOpt3 = doOpt3;
		population.resize(populationSize);
		firsttime =1;
        dist = d;

        for(int i = 0; i < populationSize; i++){
			Chromo p(noNodes , dist);
            population[i] = p;
        }

        fitness.resize(populationSize); 
        bestSoFarTour.resize(noNodes);
        for(int i = 0; i < noNodes; i++)
            bestSoFarTour[i] = i;
       iteration = 0;
        opt = new LocalSearch(dist);
    }
    
//Chromo GeneticSystem::getChromosome(int idx)
//{
 //    return population[idx];
//}
   
void GeneticSystem::setIteration(int iter)
{
        iteration = iter;
}
    
int GeneticSystem::getIteration()
{
        return iteration;
}
    
void GeneticSystem::initPopulation()
{
        for(int i = 0; i < populationSize; i++)
		{
  			population[i].shuffle();
        }

}

    void GeneticSystem::updateBestSoFarTour(){
        for(int i = 0; i < populationSize; i++)
		{
			if(computeTourLength(bestSoFarTour) > population[i].getCost() || this->firsttime ==1)
			{
				firsttime = 0;
                for(int j = 0; j < noNodes; j++)
                    bestSoFarTour[j] = population[i].getGene(j);
            }

        }
    }
std::vector<int> &GeneticSystem::getBestTour()
{
        sortPopulationAscending(population);
        return population[0].getTour();
}
    int GeneticSystem::getBestTourLength(){
		 sortPopulationAscending(population);
        return population[0].getCost();
    }
    int GeneticSystem::computeTourLength(std::vector<int> &tour){
        int len = 0;
        for(int i = 0; i < noNodes-1; i++){
            len+=(*dist)[tour[i]][tour[i+1]];
        }
		len+=(*dist)[tour[noNodes-1]][tour[0]];
        return len;
    }
  
 


	 // Selects candidate tour for crossover
int GeneticSystem::tournamentSelection()
{
	double sumFitness = 0.0;
	int idx1 = 0;
    for(int i = 0; i < populationSize; i++)
        sumFitness += fitness[i];
	double _rnd = rand();
	if(_rnd==0)
		_rnd++;
	double prob = ((double)_rnd / RAND_MAX) * sumFitness;
		
    double sum = fitness[0];
    while(sum < prob)
	{
        idx1++;
        sum+=fitness[idx1];
    }
        
	return idx1;	
        
}


	 void GeneticSystem::stepGeneration2()
	 {
		
		

		 std::vector<Chromo> newPopulation(populationSize*2);
		 for(int i = 0; i < populationSize*2; i++){
			Chromo p(noNodes , dist);
            newPopulation[i] = p;
        }

        int newPopulationSize = 0;
        sortPopulationDescending(population);
        computeFitness();

         for(int i = 0; i < populationSize; i++){
			//Chromo *Test = new Chromo(noNodes, dist);
          //  newPopulation[i] =Test;
            for(int j = 0; j < noNodes; j++)
			{
				int test = population[i].getGene(j);
                newPopulation[i].setGene(j,test );
            }
            newPopulationSize++;
        }
     

		 while (newPopulationSize < 2*populationSize) 
		 {
             // selectia 2 cromozomi pentru combinare
            //Chromo offspring1, offspring2;
			int idx1 = tournamentSelection();
			int idx2 = tournamentSelection();
         	Chromo offspring = population[idx1].CrossOver2(&population[idx2]);

            newPopulation[newPopulationSize] = offspring;
            newPopulationSize++;
        }
  
     

        mutatePopulation(newPopulation);
        sortPopulationAscending(newPopulation);
        for(int i = 0; i < populationSize; i++)
		{
          //  Chromo *p =  population[i];
			population[i] = newPopulation[i];
		//	delete p;
		}
	    sortPopulationAscending(population);
        updateBestSoFarTour();
    }



    void GeneticSystem::stepGeneration(){
		//stepGeneration2();
		//return;



		std::vector<Chromo> newPopulation(populationSize*2);
		 for(int i = 0; i < populationSize*2; i++){
			Chromo p(noNodes , dist);
            newPopulation[i] = p;
        }


        int newPopulationSize = 0;
        sortPopulationDescending(population);
        computeFitness();

        // Copiere populatie initiala in noua populatie
        for(int i = 0; i < populationSize; i++){
			//Chromo Test = new Chromo(noNodes, dist);
           // newPopulation[i] =Test;
            for(int j = 0; j < noNodes; j++)
			{
				int test = population[i].getGene(j);
                newPopulation[i].setGene(j,test );
            }
            newPopulationSize++;
        }
   //     double sumFitness = 0.0;
    //    for(int i = 0; i < populationSize; i++)
    //        sumFitness += fitness[i];

  //   Chromo *p_offspring1; 
  //   Chromo *p_offspring2;

        while(newPopulationSize < 2*populationSize)
		{
            // selectia 2 cromozomi pentru combinare
            //Chromo offspring1, offspring2;
			//p_offspring1 =new Chromo(noNodes,dist);
		//	p_offspring2= new Chromo(noNodes,dist);
       

           int idx1 =0; //tournamentSelection();
		   int idx2 =0; //tournamentSelection();
		//	printf("%d=%d,",idx1,idx2);
			while(idx1 == idx2)
			{
				idx2= tournamentSelection();
				idx1= tournamentSelection();
				
				//printf("%d=%d,",idx1,idx2);
			}
			Chromo &pfather=population[idx2];
			Chromo &pMother=population[idx1];

            Chromo p_offspring1(noNodes,dist) , p_offspring2(noNodes,dist);
			pMother.crossover(&pfather, &p_offspring1, &p_offspring2);
			// Chromo p_offspring1 = pMother.CrossOver2(&population[idx2]);
			// Chromo p_offspring2 = pfather.CrossOver2(&population[idx1]);
            newPopulation[newPopulationSize] = p_offspring1;
            newPopulationSize++;
			if(newPopulationSize >= newPopulation.size())
                break;
            newPopulation[newPopulationSize] = p_offspring2;
            newPopulationSize++;

        }


        mutatePopulation(newPopulation);
        sortPopulationAscending(newPopulation);
        for(int i = 0; i < populationSize; i++)
		{
            Chromo p =  population[i];
			
			population[i] = newPopulation[i];
		//	delete p;
		}
	

        sortPopulationAscending(population);
        updateBestSoFarTour();
    }
    void GeneticSystem::computeFitness(){
        double selectionPressure = 2;
        for(int i = 0; i < populationSize; i++){
            double fit = (double) 1*(i+1)/ (double) (populationSize-1);
            fitness[i] = fit;
        }
    }
	

    void GeneticSystem::mutatePopulation(std::vector<Chromo> &popN){
        for(int i = 0; i < popN.size(); i++)
		{
			
			double rnd =((double) rand() / (RAND_MAX)) ;
			//double prob = new java.util.Random().nextDouble()*sumFitness;
		
            if(rnd < mutationRatio)
                popN[i].mutate();
        }
    }
    void GeneticSystem::sortPopulationAscending(std::vector<Chromo> &pops){
        bool gata = false;
        while(!gata){
            gata = true;
			for(int i = 0; i < pops.size()-1; i++){
                if(pops[i].getCost() > pops[i+1].getCost()){
                    Chromo swap =pops[i];
                    pops[i] = pops[i+1];
                    pops[i+1] = swap;
                    gata = false;
                }
            }
        }
    }
    void GeneticSystem::sortPopulationDescending( std::vector<Chromo> &pops){
        bool gata = false;
        while(!gata){
            gata = true;
            for(int i = 0; i < pops.size()-1; i++){
                if(pops[i].getCost() < pops[i+1].getCost()){
                    Chromo swap=pops[i];
                    pops[i] = pops[i+1];
                    pops[i+1] = swap;
                    gata = false;
                }
            }
        }
    }
     std::vector<int> &GeneticSystem::getBestSoFarTour()
	 {
        return bestSoFarTour;
    }
    int GeneticSystem::getBestSoFarTourLength(){
        return computeTourLength(bestSoFarTour);
    }
    
	void GeneticSystem::localSearch(){
        if(doOpt2){
            for(int i = 0; i < populationSize; i++){
                opt->opt2(population[i].getTour());
            }
        }
        if(doOpt3){
            for(int i = 0; i < populationSize; i++){
                opt->opt3(population[i].getTour());
            }
        }

        updateBestSoFarTour();
    }

