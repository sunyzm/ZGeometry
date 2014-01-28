#ifndef ZUTIL_ZASSERT_H
#define ZUTIL_ZASSERT_H

#include <string>
#include <stdexcept>

namespace ZUtil {
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
} 

#endif