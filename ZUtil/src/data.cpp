#include <data.h>
#include <fstream>

using namespace std;

void vector2file( const std::string& filepath, const std::vector<double>& data, bool transpose /* = false */ )
{
	std::ofstream ofs(filepath.c_str());
	if (!transpose)
	{
		for (std::vector<double>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << ' ';
	}
	else
	{
		for (std::vector<double>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << std::endl;
	}
	ofs.close();
}

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

std::vector<double> splitStringToDouble( std::string s )
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
