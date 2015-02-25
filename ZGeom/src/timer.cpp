#include "timer.h"
#include <iostream>

using namespace std::chrono;

CStopWatch::CStopWatch()
{
    t_begin = system_clock::now();
}

void CStopWatch::startTimer()
{
    t_start = system_clock::now();
}

void CStopWatch::stopTimer()
{
    t_stop = system_clock::now();
}

void CStopWatch::stopTimer(const std::string& lead, const std::string suffix, std::ostream& os)
{
    stopTimer();
	double sec = getElapsedTime();
	os << lead << sec << suffix << std::endl;
}

double CStopWatch::getElapsedTime() const 
{
    duration<double> elapsed_seconds = t_stop - t_start;
    return elapsed_seconds.count();
}

void CStopWatch::total( const std::string& lead, std::ostream& os /*= std::cout*/ )
{
    t_end = system_clock::now();
    duration<double> elapsed_sec = t_end - t_begin;
    double sec = elapsed_sec.count();
	os << lead << sec << std::endl;
}
