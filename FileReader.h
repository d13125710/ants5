#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <array> //for std::array


#include "Tokenizer.h"
#include "TokenFinder.h"

#include "DistanceHelper.h"


using namespace std;





typedef Tokenizer<TokenFinder> Tk;
typedef std::vector<std::vector<int> > MatrixArrayType;


class FileReader
{
	struct Coords{
		double x;
		double y;
	};

private: 

	std::ifstream m_inFile; 
	std::string m_fileName;
	string sInput;
	vector<Coords> m_strVector;
 	MatrixArrayType m_theMatrix;
	IDistance *p_mDistance;

private:
	//open
	bool Openfile();
	void matrix_calculate();

public:
	bool Read(int setStartCity=0);
	double ShowSolution();
	double ReadMatrix(int x, int y) const
	{
		return  m_theMatrix[x][y];
	}
	MatrixArrayType CopyMatrix()
	{
		return m_theMatrix; 
	}
	int getsize() const
	{
	     return  (int)m_theMatrix.size();
	}
	FileReader(string &str);
	~FileReader(void); 

};

