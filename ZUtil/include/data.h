#pragma once
#include <vector>
#include <string>

void vector2file( const std::string& filepath, const std::vector<double>& data, bool transpose = false );

std::vector<int> splitStringToInt(std::string s);

std::vector<double> splitStringToDouble(std::string s);