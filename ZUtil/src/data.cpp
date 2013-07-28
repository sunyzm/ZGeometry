#include "data.h"
#include <fstream>

using namespace std;

std::vector<int> ZUtil::splitStringToInt( const std::string& s )
{
	vector<int> vd;
	size_t i,j;
	i = 0; j = s.find(',');
	for (; i != string::npos; )
	{
		if (i < j)
			vd.push_back(atoi(s.substr(i,j-i).c_str()));

		if (j == string::npos)
			i = j;
		else
		{
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;
}

std::vector<double> ZUtil::splitStringToDouble( const std::string& s )
{
	vector<double> vd;
	size_t i,j;
	i = 0; j = s.find(',');
	for (; i != string::npos; )
	{
		if (i < j)
			vd.push_back(atof(s.substr(i,j-i).c_str()));

		if (j == string::npos)
			i = j;
		else
		{
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;

}
