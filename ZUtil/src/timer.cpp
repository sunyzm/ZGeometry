#include "timer.h"
#include <iostream>

CStopWatch::CStopWatch() {
	start.QuadPart = 0;
	stop.QuadPart = 0; 
	QueryPerformanceFrequency(&frequency) ;
}

void CStopWatch::startTimer() {
	QueryPerformanceCounter(&start) ;
}

void CStopWatch::stopTimer() {
	QueryPerformanceCounter(&stop);
}

void CStopWatch::stopTimer(const std::string& lead, std::ostream& os) {
	QueryPerformanceCounter(&stop);
	double sec = getElapsedTime();
	os << lead << sec << std::endl;
}

double CStopWatch::getElapsedTime() const {
	LARGE_INTEGER time;
	time.QuadPart = stop.QuadPart - start.QuadPart;
	return LIToSecs(time) ;
}

double CStopWatch::LIToSecs(LARGE_INTEGER & L) const
{
	return (double)L.QuadPart / (double)frequency.QuadPart;
}