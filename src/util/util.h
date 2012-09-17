#pragma once
#include <string>
#include <sstream>
#include <cmath>
#include "color.h"

const double PI = 3.14159265358979323846;

inline std::string Int2String(int i)
{
	std::ostringstream ostr;
	ostr << i << std::flush;
	return ostr.str();
}

inline std::string Double2String(double f)
{
	std::ostringstream ostr;
	ostr << f << std::flush;
	return ostr.str();
}

inline double sinc(double x)
{
	if (fabs(x)<1e-8) return 1.0;
	else return std::sin(PI*x) / (PI*x);
}

class MyNote
{
public:
	int m_idx1, m_idx2;			//here id is vertex index on a particular level
	double m_score;
public:
	MyNote(int mid, double s) { m_idx1 = mid; m_score = s; }
	MyNote(int mid1, int mid2, double s) { m_idx1 = mid1; m_idx2 = mid2; m_score = s; }
	MyNote& operator = (const MyNote& note) { m_idx1 = note.m_idx1; m_idx2 = note.m_idx2; m_score = note.m_score; return(*this); }
};

class NoteCompare
{
public:
	bool operator()(const MyNote& Left, const MyNote& Right) const { return ( Left.m_score < Right.m_score ); }
};

typedef std::priority_queue<MyNote, std::vector<MyNote>, NoteCompare> NoteQueue;