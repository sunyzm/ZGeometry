#pragma once
#include <string>

class MatchTestCase
{
public:
	std::string mesh1File, mesh2File;
	std::string groundTruthFile;
	std::string presetMatchedFeatureFile;
};