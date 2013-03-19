#pragma once
#include <string>
#include <map>

class SimpleConfigLoader
{
public:
	SimpleConfigLoader(const std::string& config_file_path) { load(config_file_path); }
	void load(const std::string& config_file_path);
	void save(const std::string& config_file_path) const;
	void reload(const std::string& config_file_path);
	std::string getConfigValue(const std::string& query_string) const;
	int getConfigValueInt(const std::string& query_string) const;
	double getConfigValueDouble(const std::string& query_string) const;
private:
	std::map<std::string, std::string> m_config;
};
