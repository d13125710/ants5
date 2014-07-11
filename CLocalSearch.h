#pragma once

#include <vector>

typedef std::vector<std::vector<int> > MatrixArrayTypeInt;

class CLocalSearch
{
	MatrixArrayTypeInt  *m_Vdistance;

public:

	CLocalSearch()
	{
		m_Vdistance = 0;

	}

	CLocalSearch(MatrixArrayTypeInt  *Vdistance)
	{
		m_Vdistance = Vdistance;

	}
	void setMatrix(MatrixArrayTypeInt  *Vdistance)
	{
		if(m_Vdistance == 0)
		{
			m_Vdistance = Vdistance;
		}

	}


	void opt2(std::vector<int> &tour)
	 {
        int i;
        int a1,a2,a3,b1,b2,b3,swap;
        for(i = 1; i < tour.size()-2; i++)
		{
            a1 = (*m_Vdistance)[tour[i-1]][tour[i]];
            a2 = (*m_Vdistance)[tour[i]][tour[i+1]];
            a3 = (*m_Vdistance)[tour[i+1]][tour[i+2]];
            b1 = (*m_Vdistance)[tour[i-1]][tour[i+1]];
            b2 = (*m_Vdistance)[tour[i+1]][tour[i]];
            b3 = (*m_Vdistance)[tour[i]][tour[i+2]];
            if(a1 + a2 + a3 > b1 + b2 + b3)
			{
                swap = tour[i];
                tour[i] = tour[i+1];
                tour[i+1] = swap;
            }
		}
    }

	   void opt3(std::vector<int> &tour){
        int i,j;
        int distances[]={0,0,0,0,0,0};
        for(i = 1; i < tour.size() - 3; i++){
            // ABC
            distances[0] = (*m_Vdistance)[tour[i-1]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i+3]];
            // ACB
            distances[1] = (*m_Vdistance)[tour[i-1]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i+3]];
            // BAC
            distances[2] = (*m_Vdistance)[tour[i-1]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i+3]];
            // BCA
            distances[3] = (*m_Vdistance)[tour[i-1]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+3]];
            // CAB
            distances[4] = (*m_Vdistance)[tour[i-1]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i+3]];
            // CBA
            distances[5] = (*m_Vdistance)[tour[i-1]][tour[i+2]]+(*m_Vdistance)[tour[i+2]][tour[i+1]]+(*m_Vdistance)[tour[i+1]][tour[i]]+(*m_Vdistance)[tour[i]][tour[i+3]];
            // caut min
            int min =  std::numeric_limits<int>::max();
            int minIdx = -1;
            for(j = 0; j < 6; j++){
                if (min > distances[j]){
                    min = distances[j];
                    minIdx = j;
                }
            }
            switch(minIdx){
                case 0:     // ABC
                    break;
                case 1:     // ACB
                    swapV(tour, i+1, i+2);
                    break;
                case 2:     // BAC
                    swapV(tour, i, i+1);
                    break;
                case 3:     // BCA
                    // ABC -> BAC
                    swapV(tour, i, i+1);
                    // BAC -> BCA
                    swapV(tour, i+1, i+2);
                    break;
                case 4:     // CAB
                    // ABC -> CBA
                    swapV(tour, i, i+2);
                    // CBA -> CAB
                    swapV(tour, i+1, i+2);
                    break;
                case 5:     // CBA
                    swapV(tour, i, i+2);
                    break;

            }
            
        }

    }
    void swapV(std::vector<int> &tour , int i, int j)
	{
        int swap = tour[i];
        tour[i] = tour[j];
        tour[j] = swap;
    }

	std::vector<int> greedyPath(int noNodes)
	{
        std::vector<bool> visted(noNodes);
		std::vector<int> tour(noNodes);
		
	    int min, node;
        int i,j;

        for(i = 0; i < noNodes; i++)
            visted[i] = false;

        tour[0] = 0;
     
        visted[0] = true;

        for(i = 1; i < noNodes; i++){
            min = std::numeric_limits<int>::max();
            node = -1;
            for(j = 0; j < noNodes; j++){
                if((!visted[j])&&(j!=tour[i-1])){
                    if(min > (*m_Vdistance)[tour[i-1]][j]){
                       min = (*m_Vdistance)[tour[i-1]][j];
                       node = j;
                    }
                }
            }
            tour[i] = node;
            visted[node] = true;
        }
       return tour;
   }
	void greedyAntPath(std::vector<int> &tour)
	{
		auto noNodes = tour.size()-1;
		std::vector<bool> visted(noNodes);
			
	    int min, node;
        int i,j;

        for(i = 0; i < noNodes; i++)
            visted[i] = false;

        tour[0] = 0;
     
        visted[0] = true;

        for(i = 1; i < noNodes; i++){
            min = std::numeric_limits<int>::max();
            node = -1;
            for(j = 0; j < noNodes; j++){
                if((!visted[j])&&(j!=tour[i-1])){
                    if(min > (*m_Vdistance)[tour[i-1]][j]){
                       min = (*m_Vdistance)[tour[i-1]][j];
                       node = j;
                    }
                }
            }
            tour[i] = node;
            visted[node] = true;
        }
     
   }

	 void swap2(double *v, double *v2, int i, int j)
    /*
     * FUNCTION: auxiliary routine for sorting an integer array
     * INPUT: two arraya, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
		double tmp;

		tmp = v[i];
		v[i] = v[j];
		v[j] = tmp;
		tmp = v2[i];
		v2[i] = v2[j];
		v2[j] = tmp;
    }

    void sort2(double *v, double *v2, int left, int right)
    /*
     * FUNCTION: recursive routine (quicksort) for sorting one array; second
     * arrays does the same sequence of swaps
     * INPUT: two arrays, two indices
     * OUTPUT: none
     * (SIDE)EFFECTS: elements at position i and j of the two arrays are swapped
     */
    {
		int k, last;

		if (left >= right)
			return;
		swap2(v, v2, left, (left + right) / 2);
		last = left;
		for (k = left + 1; k <= right; k++)
			if (v[k] < v[left])
			swap2(v, v2, ++last, k);
		swap2(v, v2, left, last);
		sort2(v, v2, left, last);
		sort2(v, v2, last + 1, right);
    }




	~CLocalSearch(void)
	{
	}
};

