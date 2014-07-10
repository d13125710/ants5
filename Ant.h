#pragma once

#include <vector>


template <class RandomAccessIterator, class URNG>
void shuffle (RandomAccessIterator first, RandomAccessIterator last, URNG&& g)
{
  for (auto i=(last-first)-1; i>0; --i) {
    std::uniform_int_distribution<decltype(i)> d(0,i);
    std::swap (first[i], first[d(g)]);
  }
}



class Ant
{
	std::vector<int> m_antTour;
    std::vector<bool> m_visted;

	int noNodes;
    int tourLength;
   
public:
	int getAntTourLength() const
	{
        return tourLength;
    }
	void settourLength(int length)
	{
		tourLength=length;
	}
	bool isVisited(int node) const
	{
        return m_visted[node];
    }
	void setVisited(int node)
	{
        m_visted[node]=true;
    }
	void clearAntsMemory()
	{
		for(unsigned int i=0;i < m_visted.size();i++)
		{
				m_antTour[i]=0;
				m_visted[i]=false;
		}
	

	}
	int getTour(int node) const {
        return m_antTour[node];
    }
	std::vector<int>& getCompleteTour(){
        return m_antTour;
    }
	void setTour(int pos, int val){
        m_antTour[pos] = val;

    }
	
	Ant(int n);
	Ant(){}
	~Ant(void);
};

