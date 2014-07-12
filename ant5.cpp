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

#include "CTSPLIBFileReader.h"
#include "CAntSystem.h"
#include "CMinMaxAntSystem.h"
#include "CGeneticSystem.h"
#include "CElitistAntSystem.h"
#include "CAntColonySystem.h"
#include "BestWorstAntSystem.h"

#include "PheroMatrix.h"



//typedef std::vector<std::vector<int> > MatrixArrayTypeInt;
//typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;
//MatrixArrayTypeDouble pherom;
//MatrixArrayTypeInt ddistance;

//circle coo  (30, 5), (40, 10), (40, 20), (29, 25), (19, 25), (9, 19), (9, 9), (20, 5)
void PrintSolutionInOrder(std::vector<int> &best);
void doAnts(int num_cities , int ants ,  std::vector<std::vector<int> > *matrix , std::string file);
void doMMAnts(int num_cities , int ants , int third , std::vector<std::vector<int> > *matrix , std::string file);
void doElitest(int num_cities , int ants , int third , MatrixArrayTypeInt *matrix , std::string file);
void doGen(int num_cities , int ants , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file);
void doGen2(int ants  , int num_cities , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix);
void doACS(int ants  , int num_cities , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file);
void doBWAnts(int num_cities , int ants , std::vector<std::vector<int> > *matrix , std::string file);
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

	//test change online
	string file("berlin52");
//	string file("ulysses16");
	//string file("eil51");
///	string file("brg180");
	//string file("fl1400");
	//string file("pcb1173"); //pcb1173 pcb442 pr2392
	//string file("pcb442");
	//string file("ulysses22");
	//string file("pr2392");
//	string file("rat783");
	
	//pcb1173
	//	pcb442
	//	pr2392
	//	rat783
	//	lin318
//	string file("lin318");

	CTSPLIBFileReader fileReader(file);
	fileReader.Read();
	
	cout << "\nShortest optimized distance is = " << fileReader.ShowSolution();

	int num_cities = fileReader.getMatrixSize();


	std::vector<std::vector<int> > matrix = fileReader.CopyMatrix();

	//doBWAnts(num_cities , num_cities  , &matrix , file);

	//doAnts(num_cities , num_cities  , &matrix , file);
	for(int j=0;j<50;j++)
	doACS(num_cities  , num_cities , 0.1, 0.1, 0.9, false , false ,  &matrix , file);

	//doElitest(num_cities , num_cities , 4, &matrix , file);
	
	//AntSystem as(num_cities , num_cities);
	for(int j=0;j<50;j++)
	doMMAnts(num_cities , num_cities , 5 , &matrix , file);

	
	
	doGen( num_cities , num_cities , .9 , false , true , &matrix , file);

		

	
	system ("pause");
	return 0;
}
void doElitest(int num_cities , int ants , int third , MatrixArrayTypeInt *matrix , std::string file)
{

	cout<< "\n--ElitistAntSystem-- \n"; 
	int d =  std::numeric_limits<int>::max();

	CElitistAntSystem Ea(num_cities ,  ants , third , matrix );
	
	Ea.initPheromones();

	Ea.initAnts();
	Ea.calculateHeuristicMatrix();

	for(int i = 0; i < 5000; i++)
	{
		Ea.constructSolutions();
	//	Ea.localSearch();
		Ea.updatePheromones();
		Ea.incIteration();

		int dist= (int)Ea.getBestPathLengthSofar();
		
		int oDistance = m_Soltions[file];
		if(oDistance == dist)
		{
			cout<< " --best solution Found --" << oDistance << "at irreteration = " << i << "\n";
			break;
		}
		else if(dist < d)
		{
			d = dist;
			cout<< " best to date " << d << "\n";
		}

	}
	vector<int> citys=Ea.getBestSoFarPath();
	PrintSolutionInOrder(citys);   
	system ("pause");
}
void doACS(int num_cities  , int ants , double ro ,  double xi, double q0, bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--Ant Colony System--\n";
		CAntColonySystem acs(ants , num_cities , ro , xi , q0 , doOp1 ,doOp2,matrix);

		acs.initPheromones();
	
		acs.initAnts();
		acs.calculateHeuristicMatrix();
		cout << "Init complete\n"; 
		int d =  std::numeric_limits<int>::max();;

		for(int i = 0; i < 5000; i++)
		{
			
			acs.constructSolutions();
	 	//	acs.localSearch();
			acs.updatePheromones();
			acs.incIteration();
		
			int dist= (int)acs.getBestPathLengthSofar();
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
	vector<int> citys=acs.getBestSoFarPath();
	PrintSolutionInOrder(citys);
	
	system ("pause");



}
void doGen(int ants  , int num_cities , double third , bool doOp1 , bool doOp2 , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--GeneticSystem-- \n"; 
	int d =  std::numeric_limits<int>::max();;
	
	CGeneticSystem ga(ants , num_cities ,third, doOp1 , doOp2 , matrix);
	ga.initPopulation();

	for(int i = 0; i < 5000; i++)
	{
			ga.stepGeneration();
            ga.localSearch();
            ga.setIteration(ga.getIteration()+1);
			int dist= (int)ga.computePathLength(ga.getBestSoFarPath());
			int oDistance = m_Soltions[file];
			if(oDistance == dist)
			{
				cout<< " best solution Found " << oDistance << "at iteration = " << i << "\n";
				break;
			}
			else if(dist < d)
			{
				d = dist;
				cout<< " best to date " << d << "\n";
			}
			
	}
	vector<int> citys=ga.getBestSoFarPath();
	PrintSolutionInOrder(citys);   
	system ("pause");

}

void doMMAnts(int num_cities , int ants , int third , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "--MinMaxAntSystem-- \n"; 

	CMinMaxAntSystem mmas(num_cities , ants , third , matrix);
	int d = std::numeric_limits<int>::max();;;

	int oDistance = m_Soltions[file];
	
	mmas.initPheromones();
	mmas.calculateHeuristicMatrix();
    mmas.initAnts();
	cout << "Init complete\n"; 
	for(int i = 0; i < 50000; i++)
	{
			//cout << ".";	
			mmas.constructSolutions();
			mmas.localSearch();
			mmas.updatePheromones();
			mmas.incIteration();

			//cout << i << ",";
			//cout << "Best tour for itteration " << i << " " << as.getBestTourLength() << " \n";
			

			int dist= (int)mmas.getBestPathLengthSofar();
		
			if(oDistance == dist)
			{
				cout<< " best solution Found " << oDistance << "at iteration = " << i << "\n";
				break;
			}
			else if(dist < d)
			{ 
				d = dist;
				cout<< " best to date " << d << "\n";
			}
	}

	vector<int> citys=mmas.getBestSoFarPath();
	PrintSolutionInOrder(citys);
	system ("pause");
}
void doAnts(int num_cities , int ants , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--CAntSystem-- \n"; 

	CAntSystem as(num_cities , ants , matrix);
	int d = std::numeric_limits<int>::max();;;

	int oDistance = m_Soltions[file];

	as.initPheromones();
	as.calculateHeuristicMatrix();
	as.initAnts();
	cout << "Init complete\n"; 
	for(int i = 0; i < 5000; i++)
	{
		//cout << ".";	
		as.constructSolutions();
		// mmas.localSearch();
		as.updatePheromones();
		as.incIteration();

		//cout << i << ",";
		//cout << "Best tour for itteration " << i << " " << as.getBestTourLength() << " \n";


		int dist= (int)as.getBestPathLengthSofar();

		if(oDistance == dist)
		{
			cout<< " best solution Found " << oDistance << "at iteration = " << i << "\n";
			break;
		}
		else if(dist < d)
		{ 
			d = dist;
			cout<< " best to date " << d << "\n";
		}
	}

	vector<int> citys=as.getBestSoFarPath();
	PrintSolutionInOrder(citys);
	system ("pause");
}




void doBWAnts(int num_cities , int ants , std::vector<std::vector<int> > *matrix , std::string file)
{
	cout<< "\n--CBestWorstAntSystem-- \n"; 

	CBestWorstAntSystem as(num_cities , ants , matrix);
	int d = std::numeric_limits<int>::max();;;

	int oDistance = m_Soltions[file];

	as.initPheromones();
	as.calculateHeuristicMatrix();
	as.initAnts();
	cout << "Init complete\n"; 
	for(int i = 0; i < 5000; i++)
	{
		//cout << ".";	
		as.constructSolutions();
		//mmas.localSearch();
		as.updatePheromones();
		as.incIteration();

		//cout << i << ",";
		//cout << "Best tour for itteration " << i << " " << as.getBestTourLength() << " \n";


		int dist= (int)as.getBestPathLengthSofar();

		if(oDistance == dist)
		{
			cout<< " best solution Found " << oDistance << "at iteration = " << i << "\n";
			break;
		}
		else if(dist < d)
		{ 
			d = dist;
			cout<< " best to date " << d << "\n";
		}
	}

	vector<int> citys=as.getBestSoFarPath();
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

