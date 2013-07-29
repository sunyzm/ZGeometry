#include "timer.h"
#include <iostream>

CStopWatch::CStopWatch()
{
	mStart.QuadPart = 0;
	mStop.QuadPart = 0; 
	mEnd.QuadPart = 0;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&mBegin);
}

void CStopWatch::startTimer()
{
	QueryPerformanceCounter(&mStart);
}

void CStopWatch::stopTimer()
{
	QueryPerformanceCounter(&mStop);
}

void CStopWatch::stopTimer(const std::string& lead, std::ostream& os)
{
	QueryPerformanceCounter(&mStop);
	double sec = getElapsedTime();
	os << lead << sec << std::endl;
}

double CStopWatch::getElapsedTime() const 
{
	LARGE_INTEGER time;
	time.QuadPart = mStop.QuadPart - mStart.QuadPart;
	return LIToSecs(time) ;
}

double CStopWatch::LIToSecs(LARGE_INTEGER L) const
{
	return (double)L.QuadPart / (double)frequency.QuadPart;
}

void CStopWatch::total( const std::string& lead, std::ostream& os /*= std::cout*/ )
{
	QueryPerformanceCounter(&mEnd);
	LARGE_INTEGER time;
	time.QuadPart = mEnd.QuadPart - mBegin.QuadPart;
	double sec = LIToSecs(time);
	os << lead << sec << std::endl;
}
