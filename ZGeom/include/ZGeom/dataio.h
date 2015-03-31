#ifndef ZGEOM_DATAIO_H
#define ZGEOM_DATAIO_H
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <algorithm>

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

template<class T>
inline void vector2file(const std::string& filepath, const std::vector<T>& data, bool transpose = false)
{
	std::ofstream ofs(filepath.c_str());
	if (transpose)	// output as row vector
	{
		for (std::vector<T>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << ' ';
	}
	else	// output as column vector
	{
		for (std::vector<T>::const_iterator iter = data.begin(); iter != data.end(); ++iter)
			ofs << *iter << std::endl;
	}
	ofs.close();
}

bool fileExist(const std::string& filename);

template<typename T>
bool setHas(const std::set<T>& s, const T& query)
{
    return s.find(query) != s.end();
}

template<typename T>
bool setHasAll(const std::set<T>& s, const std::vector<T>& vQuery)
{
    for (const T& elem : vQuery) {
        if (s.find(elem) == s.end()) return false;
    }
    return true;
}

template<typename T>
std::set<T> setDiff(const std::set<T>& s1, const std::set<T>& s2)
{
    std::set<T> result;
    for (T elem : s1) {
        if (s2.find(elem) == s2.end()) result.insert(elem);
    }
    return result;
}

template<typename T>
std::set<T> setCombine(const std::set<T>& s1, const std::set<T>& s2)
{
    std::set<T> result(s1);
    for (T elem : s2) result.insert(elem);
    return result;
}

template<typename T>
bool setOverlap(const std::set<T>& s1, const std::set<T>& s2)
{
    for (const T& elem : s1)
        if (setHas<T>(s2, elem)) return true;
    return false;
}

template<class fwditer>
fwditer random_unique(fwditer begin, fwditer end, size_t num_random) 
{
    size_t left = std::distance(begin, end);
    while (num_random--) {
        fwditer r = begin;
        std::advance(r, rand() % left);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

#endif