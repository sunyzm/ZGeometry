#ifndef ZUTIL_DATA_H
#define ZUTIL_DATA_H

#include <vector>
#include <string>

namespace ZUtil
{

std::vector<int> splitStringToInt(const std::string& s);
std::vector<double> splitStringToDouble(const std::string& s);

}

#endif