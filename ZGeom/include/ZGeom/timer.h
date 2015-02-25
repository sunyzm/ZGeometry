#ifndef ZGEOM_TIMER_H
#define ZGEOM_TIMER_H
#include <string>
#include <iostream>
#include <chrono>

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
    std::chrono::time_point<std::chrono::system_clock> t_start;
    std::chrono::time_point<std::chrono::system_clock> t_stop;
    std::chrono::time_point<std::chrono::system_clock> t_begin;
    std::chrono::time_point<std::chrono::system_clock> t_end;
};


// Calls the provided work function and returns the number of seconds  
// that it takes to call that function. 
template <class Function>
inline double time_call_sec(Function&& f)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    f();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    return elapsed_seconds.count();
}

#endif