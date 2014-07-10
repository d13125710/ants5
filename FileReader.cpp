
#include "StdAfx.h"
#include "FileReader.h"




std::string::size_type sz;     // alias of size_t
 

FileReader::FileReader(string &str)
{
	m_fileName = str;
	m_inFile.open(str + ".tsp");
	p_mDistance =0;

}
double FileReader::ShowSolution()
{
	std::ifstream m_solutionFile; 
	m_solutionFile.open(m_fileName + ".opt.tour");
	bool Read = false;
	std::vector<int>  strVector;
	double distance = 0;

	if(m_solutionFile.is_open())
	{
		m_solutionFile.setf(ios::skipws);

		while(!m_solutionFile.eof() && m_solutionFile.good())
		{
			//we need to read up to the newline
			getline(m_solutionFile , sInput , '\n');
						
			if(Read)
			{
				//split up the input into tokens of whit spaces
				Tk ti(sInput.c_str() , TokenFinder(" ")); // find white spaces
						
				while( ti != Tk()) // skip the first as format in "place" , distance , distance , dont need the first
				{
					//convert the string to a double
					int position = std::stoi(*ti,&sz);
					if(position == -1)
					{
						Read = false;
						break;
					}
					strVector.push_back(position - 1); //our arrays start at 0
					++ti;
				}
				
			}
			//start of the COORDS in the TSP format , so read the coords from here....
			if(sInput == "TOUR_SECTION")
				Read = true;
		}
		m_solutionFile.close();	
		int start, end;

		for(unsigned int i = 1; i < strVector.size() ; i++)
		{
			start = strVector[i-1];
			end = strVector[i];
			distance += ReadMatrix(start , end);
		}
		//get the last city to first..
		distance += ReadMatrix(end, strVector[0]);
	
	}


	return distance;
}

bool FileReader::Read(int setStartCity)
{
	bool status = true;
	bool ReadCoords = false;

	if(m_inFile.is_open())
	{
		m_inFile.setf(ios::skipws);
		while(!m_inFile.eof() && m_inFile.good())
		{
			//we need to read up to the newline
			getline(m_inFile , sInput , '\n');
			if(ReadCoords)
			{
				//split up the input into tokens of whit spaces
				Tk ti(sInput.c_str() , TokenFinder(" ")); // find white spaces
				Coords cd;
				while( ti++ != Tk()) // skip the first as format in "place" , distance , distance , dont need the first
				{
					//convert the string to a double
					cd.x = std::stod (*ti,&sz);
					++ti; //get the next token
					cd.y = std::stod (*ti,&sz);
					//add the struct to the vector to to be stored
					m_strVector.push_back(cd);
					++ti;
				}
			}
			//start of the COORDS in the TSP format , so read the coords from here....
			std::size_t found = sInput.find("EDGE_WEIGHT_TYPE");
	
			if(found!=std::string::npos)
			{
				//find the type ..
				Tk ti1(sInput.c_str() , TokenFinder(":")); // find white spaces
				ti1++;
				if(*ti1 == "EUC_2D")
					p_mDistance =new RoundDistance();
					//round_distance
				else if(*ti1 == "CEIL_2D")
					p_mDistance = new EuclidianDistance();
					//calculateEuclidianDistance
				else if(*ti1 == "GEO")//geo_distance
					p_mDistance = new geoDistance();
				else if(*ti1 == "ATT")//att_distance
					p_mDistance = new attDistance();
				else
					p_mDistance = new EuclidianDistance();
			}
			if(sInput == "NODE_COORD_SECTION" || sInput == "EDGE_WEIGHT_SECTION")
				ReadCoords = true;
		}
		matrix_calculate();
	}
	return status;

}
void FileReader::matrix_calculate()
{
	if(p_mDistance)
	{
		size_t Vsize = m_strVector.size();
		m_theMatrix.resize(m_strVector.size());
		for(size_t i = 0 ; i< Vsize ; i++)
			m_theMatrix[i].resize(Vsize);
		
		for (int i=0; i<m_strVector.size(); i++) 
		{
			for (int j=0; j<m_strVector.size(); j++)
				m_theMatrix[i][j]  = (int) p_mDistance->calculate(m_strVector[i].x, m_strVector[i].y, m_strVector[j].x, m_strVector[j].y);
			
		}
	}
}

FileReader::~FileReader(void)
{
	m_inFile.close();
	if(p_mDistance)
		delete p_mDistance;
}

