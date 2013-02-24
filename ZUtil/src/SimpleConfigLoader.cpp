#include <SimpleConfigLoader.h>
#include <fstream>
#include <sstream>
using namespace std;

void SimpleConfigLoader::load( const std::string& config_file_path )
{
	ifstream ifs(config_file_path.c_str());
	while (ifs)
	{
		string sline;
		getline(ifs, sline);
		if (sline[0] == '#') continue;
		if (*sline.begin() == '[' && *sline.end() == ']')
		{
			string s_key = sline.substr(1, sline.size()-2);
			string s_value;
			getline(ifs, s_value);
			if (sline[0] == '#') continue;
			m_config.insert(make_pair(s_key, s_value));
		}		
	}
	ifs.close();
}

void SimpleConfigLoader::save( const std::string& config_file_path ) const
{
	ofstream ofs(config_file_path.c_str(), ios::trunc);
	for (auto iter = m_config.begin(); iter != m_config.end(); ++iter)
	{
		ofs << "[" + iter->first + "]" << endl;
		ofs << iter->second << endl;
	}
	ofs.close();
}

std::string SimpleConfigLoader::getConfigValue( const std::string& query_string ) const
{
	auto iter_result = m_config.find(query_string);
	if ( iter_result != m_config.end())
		return iter_result->second;
	else return std::string("");
}