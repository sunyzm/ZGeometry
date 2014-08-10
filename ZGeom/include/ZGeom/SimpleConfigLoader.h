#ifndef ZGEOM_SIMPLE_CONFIG_LOADER_H
#define ZGEOM_SIMPLE_CONFIG_LOADER_H
#include <string>
#include <unordered_map>

class SimpleConfigLoader
{
public:
	SimpleConfigLoader(const std::string& config_file_path) { load(config_file_path); }
	void load(const std::string& config_file_path);
	void save(const std::string& config_file_path) const;
	void reload(const std::string& config_file_path);
	std::string getConfigValue(const std::string& query_string) const;
	int			getConfigValueInt(const std::string& query_string) const;
	double		getConfigValueDouble(const std::string& query_string) const;

	bool getConfigValueString(const std::string& query_string, std::string& value_string) const;
	bool getConfigValueInt(const std::string& query_string, int& value) const;
	bool getConfigValueDouble(const std::string& query_string, double& value) const;

private:
	std::unordered_map<std::string, std::string> m_config;
};

#endif