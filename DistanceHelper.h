#pragma once
#include <math.h>
#define _USE_MATH_DEFINES

//code traken from thomas studle


const double M_PI   =        3.14159265358979323846 ; /* pi */

class IDistance
{
	public:
		IDistance(){}
		double dtrunc(double x){
			int k;

			k = (int) x;
			x = (double) k;
			return x;
		}
		virtual ~IDistance() {}
		virtual double calculate(double x1, double y1,double x2, double y2) = 0;
};


class CEuclidianDistance : public IDistance
{
	public:	
			CEuclidianDistance(){}
	double calculate(double x1, double y1,double x2, double y2) 
	{
		double x = x1 - x2;
		double y = y1 - y2;
		double dist=(x*x) + (y*y);           //calculating distance by euclidean formula
		return  sqrt(dist);    
	}
};
class CRoundDistance : public IDistance
{
	public:
		CRoundDistance() {}
	double calculate(double x1, double y1,double x2, double y2) 
	{
	   double x = x1 - x2;
	   double y = y1 - y2;
	   double dist=((x*x) + (y*y));           //calculating distance by euclidean formula
	   double round =  sqrt(dist) + 0.5;  
	   return  round;
	}

};
class CGeoDistance : public IDistance
{
	public:
		CGeoDistance(){}
	double calculate(double x1, double y1,double x2, double y2){  
		double deg, min;
		double lati, latj, longi, longj;
		double q1, q2, q3;
		double dd;
	
		deg = dtrunc(x1);
		min = x1 - deg;
		lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(x2);
		min = x2 - deg;
		latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

		deg = dtrunc(y1);
		min = y1 - deg;
		longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
		deg = dtrunc(y2);
		min = y2 - deg;
		longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

		q1 = cos(longi - longj);
		q2 = cos(lati - latj);
		q3 = cos(lati + latj);
		dd = (int) (6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		return dd;
	 }
};

class CAttDistance : public IDistance
{
public:
	CAttDistance() {}
	double calculate(double x1, double y1,double x2, double y2)
	 {
	 
		double xd = x1 - x2;
		double yd =  y1 - y2;
		double rij = sqrt((xd * xd + yd * yd) / 10.0);
		double tij = dtrunc(rij);
		int dij;

		if (tij < rij)
			dij = (int) tij + 1;
		else
			dij = (int) tij;
		return dij;
	}
};	



