#include "StdAfx.h"
#include "Ant.h"


Ant::Ant(int n)
{
	this->m_antTour.resize(n+1);
	this->m_visted.resize(n);
    tourLength = 0;
    noNodes = n;

}
Ant::~Ant(void)
{
}
