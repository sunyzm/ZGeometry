#ifndef ZGEOM_TIMER_H
#define ZGEOM_TIMER_H
#include <string>
#include <windows.h>
#include <iostream>

class CStopWatch 
{
public:
	CStopWatch();
	void startTimer();
	void stopTimer();
	void stopTimer(const std::string& lead, const std::string suffix = "", std::ostream& os = std::cout);
	void total(const std::string& lead, std::ostream& os = std::cout);
	double getElapsedTime() const;

private:
	LARGE_INTEGER mStart;
	LARGE_INTEGER mStop;
	LARGE_INTEGER mBegin;
	LARGE_INTEGER mEnd;
	LARGE_INTEGER frequency;
	double LIToSecs(LARGE_INTEGER L) const;
};


// Calls the provided work function and returns the number of milliseconds  
// that it takes to call that function. 
template <class Function> 
inline __int64 time_call(Function&& f)
{
	__int64 begin = GetTickCount64();
	f();
	return GetTickCount64() - begin;
}

template <class Function>
inline double time_call_sec(Function&& f)
{
	__int64 begin = GetTickCount64();
	f();
	return (GetTickCount64() - begin) / 1e3;
}

#endif