#include <util/output.h>
#include <fstream>


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
