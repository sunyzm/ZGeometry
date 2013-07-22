#include "timer.h"
#include <iostream>

double CStopWatch::LIToSecs( LARGE_INTEGER & L) const {
	return (double)L.QuadPart / (double)frequency.QuadPart;
}

CStopWatch::CStopWatch() {
	timer.start.QuadPart=0;
	timer.stop.QuadPart=0; 
	QueryPerformanceFrequency(&frequency) ;
}

void CStopWatch::startTimer( ) {
	QueryPerformanceCounter(&timer.start) ;
}

void CStopWatch::stopTimer( ) {
	QueryPerformanceCounter(&timer.stop);
}

void CStopWatch::stopTimer(const std::string& lead, std::ostream& os) {
	QueryPerformanceCounter(&timer.stop);
	double sec = getElapsedTime();
	os << lead << sec << std::endl;
}

double CStopWatch::getElapsedTime() const {
	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
	return LIToSecs(time) ;
}