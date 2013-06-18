#include "data.h"
#include <fstream>

using namespace std;

namespace ZUtil {

std::vector<int> splitStringToInt( std::string s )
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

std::vector<double> ZUtil::splitStringToDouble( std::string s )
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

bool fileExist( const std::string& filename )
{
	std::ifstream ifs(filename.c_str());
	bool retval = true;
	if (!ifs) retval = false;
	ifs.close();
	return retval;
}

}


