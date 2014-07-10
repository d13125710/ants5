// ant5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <vector>
#include <iostream>
#include <time.h>
#include <math.h>
#include <random>
#include <string>
#include <map>

#include "FileReader.h"
#include "AntSystem.h"
#include "MinMaxAntSystem.h"
#include "GeneticSystem.h"
#include "ElitistAntSystem.h"
#include "AntCSystem.h"



//typedef std::vector<std::vector<int> > MatrixArrayTypeInt;
//typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;
//MatrixArrayTypeDouble pherom;
//MatrixArrayTypeInt ddistance;

//circle coo  (30, 5), (40, 10), (40, 20), (29, 25), (19, 25), (9, 19), (9, 9), (20, 5)
void PrintSolutionInOrder(std::vector<int> &best);
void doAnts(int num_cities , int ants , int third , std::vector<std::vector<int> > *matrix , std::string file);
void doElitest(int num_cities , int ants , int third , MatrixArrayTypeInt *matrix , std::string file);
void doGen(int num_cities , int ants , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file);
void doGen2(int ants  , int num_cities , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix);
void doACS(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file);

std::map<std::string,unsigned int> m_Soltions;


typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;
typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

void populateSolutionMap(std::string &file  , unsigned int distance)
{
	m_Soltions[file] = distance;
}

void createSolutionsMap()
{
	populateSolutionMap(string("berlin52") , 7542);
	populateSolutionMap(string("eil51") , 426);
	populateSolutionMap(string("pr2392") , 378032);
	populateSolutionMap(string("ulysses16") , 6859);
	populateSolutionMap(string("fl1400") , 20127);
	populateSolutionMap(string("rat783") , 8806);
	populateSolutionMap(string("pr2392") , 378032);

}



int _tmain(int argc, _TCHAR* argv[])
{
	std::srand((unsigned )time(NULL));
    srand( (unsigned)time( NULL ) );
	createSolutionsMap();


	string file("berlin52");
	//string file("ulysses16");
	//string file("eil51");
	//string file("brg180");
	//string file("fl1400");
	//string file("pcb1173"); //pcb1173 pcb442 pr2392
	//string file("pcb442");
	//string file("ulysses22");
	//string file("pr2392");
	//string file("rat783");
	
	//pcb1173
	//	pcb442
	//	pr2392
	//	rat783
	//	lin318
//	string file("lin318");

	FileReader fileReader(file);
	fileReader.Read();
	
	cout << "\nShorest optimized distance is = " << fileReader.ShowSolution();

	int num_cities = fileReader.getsize();


	std::vector<std::vector<int> > matrix = fileReader.CopyMatrix();


//	doElitest(num_cities , 20 , 4, &matrix , file);
	
	//AntSystem as(num_cities , num_cities);
	doAnts(num_cities , 100 , 10 , &matrix , file);

	// doGen2(10 , num_cities , 0.30 , false , false , matrix);

	
	// doGen(318 , num_cities , .9 , false , true , matrix , file);

		
	//doACS(num_cities  , num_cities , 0.1, 0.1, 0.9, false , false ,  matrix , file);
	

	return 0;
}
void doElitest(int num_cities , int ants , int third , MatrixArrayTypeInt *matrix , std::string file)
{

	cout<< "\n--ElitistAntSystem-- \n"; 
	int d =  std::numeric_limits<int>::max();

	ElitistAntSystem Ea(num_cities ,  ants , third , matrix );
	
		Ea.initPheromones();

	Ea.initAnts();
	Ea.calculateHeuristicMatrix();

	for(int i = 0; i < 1000; i++)
	{
		Ea.constructSolutions();
		//Ea.localSearch();
		Ea.updatePheromones();
		Ea.incIteration();

		int dist= (int)Ea.getBestTourlengthSofar();
		
		int oDistance = m_Soltions[file];
		if(oDistance == dist)
		{
			cout<< " best solution Found " << oDistance << "at irreteration = " << i << "\n";
			break;
		}
		else if(dist < d)
		{
			d = dist;
			cout<< " best to date " << d << "\n";
		}

	}
	vector<int> citys=Ea.getBestSoFarTour();
	PrintSolutionInOrder(citys);   
	system ("pause");


}


void doACS(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--Ant Colony System--\n";
		AntCSystem acs(ants , num_cities , ro , xi , q0 , doOp1 ,doOp2,matrix);

		acs.initPheromones();
	
		acs.initAnts();
		acs.calculateHeuristicMatrix();
		cout << "Init complete\n"; 
		int d =  std::numeric_limits<int>::max();;

		for(int i = 0; i < 5000; i++)
		{
			
			acs.constructSolutions();
	 		acs.localSearch();
			acs.updatePheromones();
			acs.incIteration();
		
			int dist= (int)acs.getBestTourlengthSofar();
			int oDistance = m_Soltions[file];

			if(oDistance == dist)
			{
				cout<< " best solution Found "  << oDistance << "at irreteration = " << i << "\n";
			
				break;
			}
			else if(dist < d)
			{
				d = dist;
				cout<< " best to date " << d << "\n";
			}
			
			
		}
	vector<int> citys=acs.getBestSoFarTour();
	PrintSolutionInOrder(citys);
	
	system ("pause");



}
void doGen(int ants  , int num_cities , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--GeneticSystem-- \n"; 
	int d =  std::numeric_limits<int>::max();;
	
	GeneticSystem ga(ants , num_cities ,third, doOp1 , doOp2 , matrix);
	ga.initPopulation();

	for(int i = 0; i < 1000; i++)
	{
			ga.stepGeneration();
            ga.localSearch();
            ga.setIteration(ga.getIteration()+1);
			int dist= (int)ga.computeTourLength(ga.getBestSoFarTour());
			int oDistance = m_Soltions[file];
			if(oDistance == dist)
			{
				cout<< " best solution Found " << oDistance << "at irreteration = " << i << "\n";
				break;
			}
			else if(dist < d)
			{
				d = dist;
				cout<< " best to date " << d << "\n";
			}
			
	}
	vector<int> citys=ga.getBestSoFarTour();
	PrintSolutionInOrder(citys);   
	system ("pause");

}

void doAnts(int num_cities , int ants , int third , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "--MinMaxAntSystem-- \n"; 

	MinMaxAntSystem mmas(num_cities , ants , third , matrix);
	int d = std::numeric_limits<int>::max();;;

	int oDistance = m_Soltions[file];
	
	mmas.initPheromones();
	mmas.calculateHeuristicMatrix();
    mmas.initAnts();
	cout << "Init complete\n"; 
	for(int i = 0; i < 5000; i++)
	{
			//cout << ".";	
			mmas.constructSolutions();
			// mmas.localSearch();
			mmas.updatePheromones();
			mmas.incIteration();

			//cout << i << ",";
			//cout << "Best tour for itteration " << i << " " << as.getBestTourLength() << " \n";
			

			int dist= (int)mmas.getBestTourlengthSofar();
		
			if(oDistance == dist)
			{
				cout<< " best solution Found " << oDistance << "at irreteration = " << i << "\n";
				break;
			}
			else if(dist < d)
			{
				d = dist;
				cout<< " best to date " << d << "\n";
			}
	}

	vector<int> citys=mmas.getBestSoFarTour();
	PrintSolutionInOrder(citys);
	system ("pause");
}

//find the number 0
//if its the beg just print out the tour
//if its the end print tour backwards

//else find the pos
//print from that pos to the end
//print from that pos to the start

void PrintSolutionInOrder(std::vector<int> &best)
{
	vector<int>::iterator it , itend , itbeg , itpos; 
	itpos=find(best.begin(),best.end(),0);
	
	if(itpos == best.begin())
	{
		//print the string
		it=itpos;
		for (it; it != best.end(); it++) 
		{
			printf("%i|",(*it)+1);
		}

	}
	else if(itpos == best.end())
	{
		//print the string backwards
	}
	else
	{
		it = itpos;
		for (it; it != best.end(); it++) 
		{
			printf("%i|",(*it)+1);
		}
		it = best.begin();
		for (it; it != itpos; it++) 
		{
			printf("%i|",(*it)+1);
		}
	}
	printf("\n");


	
}

