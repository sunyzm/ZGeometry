#include "dataio.h"
#include <fstream>
#include <filesystem>

using namespace std;

std::vector<int> splitStringToInt( const std::string& s )
{
	vector<int> vd;
	size_t i,j;
	i = 0; j = s.find(',');
	while (i != string::npos)
	{
		if (i < j)
			vd.push_back(atoi(s.substr(i,j-i).c_str()));

		if (j == string::npos) i = j;
		else {
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;
}

std::vector<double> splitStringToDouble( const std::string& s )
{
	vector<double> vd;
	size_t i,j;
	i = 0; j = s.find(',');
    while (i != string::npos)
	{
		if (i < j)
			vd.push_back(atof(s.substr(i,j-i).c_str()));

		if (j == string::npos) i = j;
		else {
			i = j + 1;
			j = s.find(',',j+1);
		}
	}
	return vd;
}

bool fileExist(const std::string& filename)
{
    std::tr2::sys::path myfile(filename);
    return std::tr2::sys::exists(myfile);
}
