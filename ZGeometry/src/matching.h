#pragma once
#include <map>
#include <vector>
#include <utility>
#include <string>
#include <functional>
#include <iostream>

class MatchPair
{
public:
	int		m_idx1;
	int		m_idx2;
	double	m_tl;		// lower time
	double	m_tu;		// upper time
	int		m_tn;		// number of timescales
	double  m_score;	
	int		m_note;		

public:
	MatchPair() { m_idx1 = -1; m_idx2 = -1; m_score = m_tl = m_tu = 0;m_tn = m_note = 0; }
	MatchPair(int i1, int i2) : m_idx1(i1), m_idx2(i2) { m_score = m_tl = m_tu = 0; m_tn = m_note = 0; }
	MatchPair(int i1, int i2, double score) : m_idx1(i1), m_idx2(i2), m_score(score) { m_tl = m_tu = 0.0; m_tn = m_note = 0;}
	MatchPair(int i1, int i2, double tl, int tn, double score = 0) :m_idx1(i1), m_idx2(i2), m_tl(tl), m_tn(tn), m_score(score) { m_tu = m_tl * pow(2., tn-1); m_note = 0;}
	operator std::pair<int, int> () const { return std::make_pair(m_idx1, m_idx2); }

	friend bool operator == (const MatchPair& mp1, const MatchPair& mp2) { return (mp1.m_idx1 == mp2.m_idx1 && mp1.m_idx2 == mp2.m_idx2); }
	friend bool operator < (const MatchPair& mp1, const MatchPair& mp2)
	{
		return mp1.m_idx1 < mp2.m_idx1 || (mp1.m_idx1 == mp2.m_idx1 && mp1.m_idx2 < mp2.m_idx2);
	}
	friend bool operator > (const MatchPair& mp1, const MatchPair& mp2)
	{
		return mp2 < mp1;
	}

	static std::vector<std::pair<int, int> > ToPairVector(const std::vector<MatchPair>& vmp);
};

const std::function<bool(const MatchPair& left, const MatchPair& right)> PairScoreCompare = [](const MatchPair& left, const MatchPair& right) { return left.m_score > right.m_score; };
/*
class PairScoreCompare
{
public:
	bool operator()(const MatchPair& Left, const MatchPair& Right) const { return ( Left.m_score > Right.m_score ); }
};
*/

class MatchResult
{
public:
	std::map<int, int> mMatchedPairs;
	void read(const std::string& file);
	void write(const std::string& file) const;
};

class MatchEvaluation
{
public:
	int mMatchedCount;
	int mTotalToMatch;
	int mMatch0Count;
	int mMatch1Count;
	int mMatch2Count;
	int mLargeErrorCount;
	double mAvgErrInEdgeLength;	// as multiples of average edge length
	double mAvgRelativeError;	// average relative errors tested from selected pairs

	void printStats( bool withGroundTruth, bool withRelativeError, std::ostream& ostr = std::cout ) const;
	void setTotal(int n) { mTotalToMatch = n; }
};