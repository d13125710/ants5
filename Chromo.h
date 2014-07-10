#pragma once

#include <vector>
#include <random>
#include <time.h>
#include <math.h>


template <class RandomAccessIterator, class URNG>
void shuffleVector (RandomAccessIterator first, RandomAccessIterator last, URNG&& g)
{
  for (auto i=(last-first)-1; i>0; --i) {
    std::uniform_int_distribution<decltype(i)> d(0,i);
    std::swap (first[i], first[d(g)]);
  }
}


typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class Chromo
{



public:
	//cctor
	 Chromo( ){ }

	 Chromo( const Chromo& other )
	 {
		 if(this != &other)
		 {
			 m_Vgenes = other.m_Vgenes;
			 m_noGenes = other.m_noGenes;
			 m_dist = other.m_dist;
			 fitness = other.fitness;
		 }


	 }


	std::vector<int> m_Vgenes;
	int m_noGenes;
	MatrixArrayTypeInt *m_dist;
	double fitness;

	std::vector<int>& getTour();
	int getNoGenes();
	int getGene(int idx);
	void setGene(int idx, int val);
	void mutate();
	void shuffle();
	void crossover(Chromo  *father,  Chromo  *offspring1,  Chromo  *offspring2);
	Chromo CrossOver2(Chromo *pParent1);

	int getCost();
   
   
	 Chromo(int noGenes, MatrixArrayTypeInt *distance);
	 ~Chromo(void);
	 double getFitness();
	 int computeTourLength();
		int getDist(int indx1,int indx2) const
		{
			return (*m_dist)[indx1][indx2];
		}

		int getRandomNumber(int low, int high)
		{
			int temp = (std::rand() % ((high + 1) - low)) + low;
			return temp;
		}
	 
	 int getRandomNumber(int low, int high, int except)
	 {
		//Generate random value between low and high (inclusive), except 'except'.
		int temp;
		bool done = false;

		while(!done)
		{
			temp = (std::rand() % ((high + 1) - low)) + low;
			if(temp != except)
			{
				done = true;
			}
		}
		return temp;
	 }

	  bool containsGene(int gene) const 
	  {
      
		for(int i =0; i< m_Vgenes.size(); i++){
			if(m_Vgenes[i]==gene)
				return true;
		}
		return false;
	  }

};

