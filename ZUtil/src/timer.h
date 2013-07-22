#ifndef ZUTIL_TIMER_H
#define ZUTIL_TIMER_H

#include <string>
#include <windows.h>
#include <iostream>

typedef struct 
{
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} StopWatch;

class CStopWatch 
{
public:
	CStopWatch();
	void startTimer();
	void stopTimer();
	void stopTimer(const std::string& lead, std::ostream& os = std::cout);
	double getElapsedTime() const;

private:
	StopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L) const;
};

// Calls the provided work function and returns the number of milliseconds  
// that it takes to call that function. 
template <class Function>
__int64 time_call(Function&& f)
{
	__int64 begin = GetTickCount();
	f();
	return GetTickCount() - begin;
}

#endif