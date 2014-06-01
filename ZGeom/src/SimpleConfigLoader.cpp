#include "SimpleConfigLoader.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include "dataio.h"

using namespace std;

void SimpleConfigLoader::load( const std::string& config_file_path )
{
	if (!fileExist(config_file_path))
		throw std::runtime_error("Cannot found configuration files!");

	ifstream ifs(config_file_path.c_str());
	while (ifs)
	{
		string sline;
		getline(ifs, sline);	
		if (sline.empty() || sline[0] == '#' || sline[0] == ';') continue;
		sline = sline.substr(0, sline.find('#'));
		
		size_t sep = sline.find('=');
		if (sep == string::npos) continue;	// no '=' is found, the line is skipped
		string s_key = sline.substr(0, sep), s_value = sline.substr(sep+1);

		string::iterator iter = s_key.begin();
		while (iter != s_key.end())
		{
			if (*iter == ' ' || *iter == '\t') iter = s_key.erase(iter);
			else ++iter;
		}
		iter = s_value.begin();
		while (iter != s_value.end())
		{
			if (*iter == ' ' || *iter == '\t') iter = s_value.erase(iter);
			else ++iter;
		}
		m_config.insert(make_pair(s_key, s_value));
	}
	ifs.close();
}

void SimpleConfigLoader::save( const std::string& config_file_path ) const
{
	ofstream ofs(config_file_path.c_str(), ios::trunc);
	for (auto iter = m_config.begin(); iter != m_config.end(); ++iter)
	{
		ofs << iter->first << '=' << iter->second << endl;
	}
	ofs.close();
}

void SimpleConfigLoader::reload( const std::string& config_file_path )
{
	m_config.clear();
	load(config_file_path);
}

std::string SimpleConfigLoader::getConfigValue( const std::string& query_string ) const
{
	auto iter_result = m_config.find(query_string);
	if ( iter_result != m_config.end())
		return iter_result->second;
	else return std::string("");
}

int SimpleConfigLoader::getConfigValueInt( const std::string& query_string ) const
{
	return atoi(getConfigValue(query_string).c_str());	
}

bool SimpleConfigLoader::getConfigValueInt( const std::string& query_string, int& value ) const
{
	std::string strValue;
	if (getConfigValueString(query_string, strValue)) {
		value = atoi(strValue.c_str());
		return true;
	} else {
		return false;
	}
}

double SimpleConfigLoader::getConfigValueDouble( const std::string& query_string ) const
{
	return atof(getConfigValue(query_string).c_str());
}

bool SimpleConfigLoader::getConfigValueDouble( const std::string& query_string, double& value ) const
{
	std::string strValue;
	if (getConfigValueString(query_string, strValue)) {
		value = atof(strValue.c_str());
		return true;
	} else {
		return false;
	}
}

bool SimpleConfigLoader::getConfigValueString( const std::string& query_string, std::string& value_string ) const
{
	auto iter_result = m_config.find(query_string);
	if (iter_result != m_config.end()) {
		value_string = iter_result->second;
		return true;
	}
	else return false;
}

