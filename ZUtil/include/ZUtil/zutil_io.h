#ifndef ZUTIL_IO_H
#define ZUTIL_IO_H

#include <string>
#include <fstream>
#include <vector>

namespace ZUtil
{
	template<class T>
	inline void vector2file( const std::string& filepath, const std::vector<T>& data, bool transpose = false )
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

	inline bool fileExist( const std::string& filename )
	{
		std::ifstream ifs(filename.c_str());
		bool retval = true;
		if (!ifs) retval = false;
		ifs.close();
		return retval;
	}
}

#endif