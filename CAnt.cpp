#include "StdAfx.h"
#include "CAnt.h"


CAnt::CAnt(int n)
{
	this->m_antTour.resize(n+1);
	this->m_visted.resize(n);
    tourLength = 0;
    noNodes = n;

}
CAnt::~CAnt(void)
{
}
