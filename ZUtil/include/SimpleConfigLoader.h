#pragma once
#include <string>
#include <map>

class SimpleConfigLoader
{
public:
	void load(const std::string& config_file_path);
	void save(const std::string& config_file_path) const;
	std::string getConfigValue(const std::string& query_string) const;
private:
	std::map<std::string, std::string> m_config;
};
