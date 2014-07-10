#include "StdAfx.h"
#include "Chromo.h"


#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector

 


Chromo::Chromo(int noGenes, MatrixArrayTypeInt *distance)
{
	this->m_noGenes= noGenes;
	m_dist = distance;

	m_Vgenes.resize(noGenes);
	fitness=0;

}

Chromo::~Chromo(void)
{

}


std::vector<int> &Chromo::getTour()
{
	return m_Vgenes;
}
int Chromo::getNoGenes()
{
    return m_noGenes;
}
int Chromo::getGene(int idx)
{
    return m_Vgenes[idx];
}
void Chromo::setGene(int idx, int val)
{
    m_Vgenes[idx] = val;
}
void Chromo::mutate()
{
    
    int idx1 = getRandomNumber(0, m_noGenes-1);
	int idx2  = getRandomNumber(0,m_noGenes-1, idx1);
	
    int swap = m_Vgenes[idx1];
    m_Vgenes[idx1] = m_Vgenes[idx2];
    m_Vgenes[idx2] = swap;
 
}
void Chromo::shuffle(){
  
    for(int i = 0; i < m_noGenes; i++)
	{
        m_Vgenes[i] = i;
    }

  std::random_shuffle( m_Vgenes.begin() , m_Vgenes.end() );

}
 Chromo Chromo::CrossOver2(Chromo *pParent1)
 {
	// Chromo *pParent1 = this;

	Chromo offspring (this->m_noGenes, this->m_dist);
	//std::vector<int> child(m_noGenes+1);
	
	for (int i = 0; i < 	offspring.m_Vgenes.size(); i++) {
			offspring.m_Vgenes[i]=-1;
	}

		
        // Get start and end sub tour positions for parent1's tour
        int startPos = (int) getRandomNumber(0 , m_noGenes);
        int endPos = (int)  getRandomNumber(0 , m_noGenes , startPos);

        // Loop and add the sub tour from parent1 to our child
		for (int i = 0; i < 	offspring.m_noGenes; i++) {
            // If our start position is less than the end position
            if (startPos < endPos && i > startPos && i < endPos)
			{
				offspring.m_Vgenes[i] = this->m_Vgenes[i];
              //  child.setCity(i, parent1.getCity(i));
            } // If our start position is larger
            else if (startPos > endPos) 
			{
                if (!(i < startPos && i > endPos)) {
                   	offspring.m_Vgenes[i] = this->m_Vgenes[i];
                }
            }
        }

        // Loop through parent2's city tour
        for (int i = 0; i < pParent1->m_noGenes; i++) {
            // If child doesn't have the city add it
			if (!offspring.containsGene(pParent1->m_Vgenes[i])) 
			{
                // Loop to find a spare position in the child's tour
				for (int ii = 0; ii < 	offspring.m_noGenes; ii++) 
				{
                    // Spare position found, add city
                    if (offspring.m_Vgenes[ii] == -1) 
					{
                        offspring.m_Vgenes[ii]=pParent1->m_Vgenes[i];
                        break;
                    }
                }
            }
        }
        return offspring;


 }


void Chromo::crossover(Chromo  *father,  Chromo *offspring1,  Chromo  *offspring2)
{
    offspring1->m_Vgenes[0] = m_Vgenes[0];
    offspring2->m_Vgenes[0] = father->m_Vgenes[0];
	int max = this->m_noGenes-1;
        
	std::vector<bool> visited1(m_noGenes), visited2(m_noGenes);

	for(int i = 0; i < m_noGenes; i++){
        visited1[i] = false;
        visited2[i] = false;
    }

    visited1[offspring1->m_Vgenes[0]] = true;
    visited2[offspring2->m_Vgenes[0]] = true;

    for(int i = 1; i < m_noGenes; i++){
        int prevNode = offspring1->m_Vgenes[i-1];
        int node1 = 0, node2 = 0;
        for(int j = 0; j < m_noGenes-1; j++){
            if(m_Vgenes[j] == prevNode)
                node1 = m_Vgenes[j+1];
            if(father->m_Vgenes[j] == prevNode)
                node2 = father->m_Vgenes[j+1];

        }
        if((!visited1[node1]) && (!visited1[node2])){
            if(getDist(prevNode,node1) < getDist(prevNode, node2)){
                offspring1->m_Vgenes[i] = node1;
                visited1[node1] = true;
            }
            else{
                offspring1->m_Vgenes[i] = node2;
                visited1[node2] = true;
            }
            continue;
        }
        if((!visited1[node1]) && visited1[node2]){
            offspring1->m_Vgenes[i] = node1;
            visited1[node1] = true;
            continue;
        }
        if((!visited1[node2]) && visited1[node1]){
            offspring1->m_Vgenes[i] = node2;
            visited1[node2] = true;
            continue;
        }
        if(visited1[node1] && visited1[node2]){
      		int node = getRandomNumber(0, max);
	        while(visited1[node]){
				
				node = getRandomNumber(0,max , node);
	      }
            visited1[node] = true;
            offspring1->m_Vgenes[i] = node;
            continue;
        }
    }
    //===================================================
    for(int i = 1; i < m_noGenes; i++){
        int prevNode = offspring2->m_Vgenes[i-1];
        int node1 = 0, node2 = 0;
        for(int j = 0; j < m_noGenes-1; j++){
            if(m_Vgenes[j] == prevNode)
                node1 = m_Vgenes[j+1];
            if(father->m_Vgenes[j] == prevNode)
                node2 = father->m_Vgenes[j+1];

        }
        if((!visited2[node1]) && (!visited2[node2])){
            if(getDist(prevNode,node1) < getDist(prevNode,node2)){
                offspring2->m_Vgenes[i] = node1;
                visited2[node1] = true;
            }
            else{
                offspring2->m_Vgenes[i] = node2;
                visited2[node2] = true;
            }
            continue;
        }
        if((!visited2[node1]) && visited2[node2]){
            offspring2->m_Vgenes[i] = node1;
            visited2[node1] = true;
            continue;
        }
        if((!visited2[node2]) && visited2[node1]){
            offspring2->m_Vgenes[i] = node2;
            visited2[node2] = true;
            continue;
        }
        if(visited2[node1] && visited2[node2]){
      	int node = getRandomNumber(0, max);
				while(visited2[node]){
				node = getRandomNumber(0,max , node);
	       }
            visited2[node] = true;
            offspring2->m_Vgenes[i] = node;
            continue;
        }

    }
  
}

int Chromo::getCost()
{
    if(fitness)
		return (int)fitness;
	
	int sum = 0;

    for(int i = 0; i < m_noGenes-1; i++){
        sum+=(*m_dist)[m_Vgenes[i]][m_Vgenes[i+1]];
    }
	sum += (*m_dist)[m_Vgenes[m_Vgenes.size()-1]][m_Vgenes[0]];
	if(this->fitness=0)
		fitness =sum;

    return sum;
}





