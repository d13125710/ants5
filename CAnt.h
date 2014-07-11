#pragma once

#include <vector>


class CAnt
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
	void setAntTourLength(int length)
	{
		tourLength=length;
	}
	bool isCityVisited(int city) const
	{
        return m_visted[city];
    }
	void setCityVisited(int city)
	{
        m_visted[city]=true;
    }
	void ReinitiseAnt()
	{
		for(unsigned int i=0;i < m_visted.size();i++)
		{
				m_antTour[i]=0;
				m_visted[i]=false;
		}
	

	}
	int getCity(int index) const {
        return m_antTour[index];
    }
	std::vector<int>& getAntsCityTour(){
        return m_antTour;
    }
	void setAntCity(int index, int city){
        m_antTour[index] = city;

    }
	
	CAnt(int noOfCitys);
	CAnt(){}
	~CAnt(void);
};

