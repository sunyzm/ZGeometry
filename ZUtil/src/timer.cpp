//#include <windows.h>
#include <timer.h>


double CStopWatch::LIToSecs( LARGE_INTEGER & L) const {
	return (double)L.QuadPart / (double)frequency.QuadPart;
}

CStopWatch::CStopWatch(){
	timer.start.QuadPart=0;
	timer.stop.QuadPart=0; 
	QueryPerformanceFrequency(&frequency) ;
}

void CStopWatch::startTimer( ) {
	QueryPerformanceCounter(&timer.start) ;
}

void CStopWatch::stopTimer( ) {
	QueryPerformanceCounter(&timer.stop) ;
}

double CStopWatch::getElapsedTime() const {
	LARGE_INTEGER time;
	time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
	return LIToSecs(time) ;
}