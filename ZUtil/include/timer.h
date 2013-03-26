#pragma once
#include <windows.h>
#include <string>

typedef struct 
{
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} StopWatch;

class CStopWatch 
{
private:
	StopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L) const;
public:
	CStopWatch();
	void startTimer( );
	void stopTimer( );
	double getElapsedTime() const;
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