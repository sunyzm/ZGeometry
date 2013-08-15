#ifndef ZUTIL_DATA_H
#define ZUTIL_DATA_H

#include <vector>
#include <string>
#include <sstream>

namespace ZUtil
{

std::vector<int> splitStringToInt(const std::string& s);
std::vector<double> splitStringToDouble(const std::string& s);

inline std::string Int2String(int i)
{
    std::ostringstream ostr;
    ostr << i << std::flush;
    return ostr.str();
}

inline std::string Double2String(double f)
{
    std::ostringstream ostr;
    ostr << f << std::flush;
    return ostr.str();
}

} //end of namespace ZUtil

#endif