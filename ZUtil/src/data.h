#ifndef ZUTIL_DATA_H

#include <vector>
#include <string>

template<class T>
void vector2file( const std::string& filepath, const std::vector<T>& data, bool transpose = false )
{
	std::ofstream ofs(filepath.c_str());
	if (transpose)	// output as row vector
	{
		for (std::vector<T>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << ' ';
	}
	else	// output as column vector
	{
		for (std::vector<T>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << std::endl;
	}
	ofs.close();
}

std::vector<int> splitStringToInt(std::string s);

std::vector<double> splitStringToDouble(std::string s);

#endif