#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand





typedef std::vector<std::vector<int> > MatrixArrayTypeInt;


 static int myrandom2 (int i) { return std::rand()%i;}

class Chromo1
{
	std::vector<int> tour;
    int m_noGenes;
	MatrixArrayTypeInt *m_dist;
	
	double fitness ;
    int distance ;
	

public:
	int GetSize()
	{
		return m_noGenes;
	}
	 void generateIndividual() 
	 {
		//  std::srand ( unsigned ( std::time(0) ) );
        // srand((unsigned int)time(0));
		 for(int i = 0; i < m_noGenes; i++)
			   tour[i] = i;
  
		//std::mt19937   mt_rand(time(NULL));
		
		std::random_shuffle(tour.begin(), tour.end() , myrandom2);
    }
   // Gets the tours fitness
    double getFitness() 
	{
        if (fitness == 0) 
            fitness = 1/(double)computeTourLength();
        return fitness;
    }
	int computeTourLength()
	{
		int  start, end;
		double d = 0;
		for (unsigned int i=1; i< tour.size(); i++)
		{
			start =tour[i-1];
			end = tour[i];
			d += (*m_dist)[start][end];	
		}

		d += (*m_dist)[tour[tour.size()-1]][tour[0]];
		return (int)d;  
	  
	}

   void addCity(int pos , int city) {
		 tour[pos] = city;
     
    }
    // Get a city
      int getCity(int index){
        return tour[index];
    }
    
    // Get the number of destination cities
      size_t numberOfCities(){
        return tour.size();
    }
	   // Check if the tour contains a city
    bool containsCity(int city){
      
		for(int i =0; i< tour.size(); i++)
		{
			if(tour[i]==city)
				return true;
		}
		return false;
		//return tour.contains(city);
    }



	Chromo1(int noGenes, MatrixArrayTypeInt *distance);
	~Chromo1(void);
};

