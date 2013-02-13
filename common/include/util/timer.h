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

