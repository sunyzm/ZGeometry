#ifndef ZUTIL_ZASSERT_H
#define ZUTIL_ZASSERT_H

#include <string>
#include <stdexcept>

namespace ZUtil{
	inline void logic_assert(bool condition, const char* errorMsg = "logic_error") 
	{
		if (!condition) {
			throw std::logic_error(errorMsg);
		}
	}

	inline void runtime_assert(bool condition, const char* errorMsg = "runtime_error")
	{
		if (!condition) {
			throw std::runtime_error(errorMsg);
		}
	}
} 

#endif