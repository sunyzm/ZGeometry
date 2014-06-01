#ifndef ZGEOM_ZASSERT_H
#define ZGEOM_ZASSERT_H
#include <string>
#include <stdexcept>

namespace ZGeom {

inline void logic_assert(bool condition, const std::string& errorMsg = "logic_error")
{
	if (!condition) {
		throw std::logic_error(errorMsg.c_str());
	}
}

inline void runtime_assert(bool condition, const std::string& errorMsg = "runtime_error")
{
	if (!condition) {
		throw std::runtime_error(errorMsg.c_str());
	}
}

}// end of namespace

#endif