#include "matching.h"
#include <fstream>
#include <ZUtil/ZUtil.h>

std::vector<std::pair<int, int> > MatchPair::ToPairVector(const std::vector<MatchPair>& vmp)
{
	std::vector<std::pair<int, int> > vp;
	for (const MatchPair& mp : vmp) vp.push_back((std::pair<int, int>)mp);
	return vp;
}

void MatchResult::write( const std::string& file ) const
{
	std::ofstream ofs(file.c_str());
	ofs << mMatchedPairs.size() << std::endl;
	for (auto p : mMatchedPairs) {
		ofs << p.first << ' ' << p.second << std::endl;
	}
	ofs.close();
}

void MatchResult::read( const std::string& file )
{
	std::ifstream ifs(file.c_str());
	ZUtil::runtime_assert(ifs, "File " + file + " not exist!");

	int size, idx1, idx2;
	ifs >> size;
	for (int i = 0; i < size; ++i) {
		ifs >> idx1 >> idx2;
		mMatchedPairs.insert(std::make_pair(idx1, idx2));
	}
	ifs.close();
}

void MatchEvaluation::printStats( bool withGroundTruth, bool withRelativeError, std::ostream& ostr /*= std::cout */ ) const
{
	ostr << "** Match Evaluation Results **" << std::endl;
	ostr << "Matched count: " << mMatchedCount << "/" << mTotalToMatch << "(" << double(mMatchedCount)/mTotalToMatch << ")" << std::endl;

	if (withGroundTruth) {
		ostr << "Precise match count: " << mMatch0Count << "(" << double(mMatch0Count)/mTotalToMatch << ")" << std::endl;
		ostr << "1-ring match count: " << mMatch1Count << "(" << double(mMatch1Count)/mTotalToMatch << ")" << std::endl;
		ostr << "2-ring match count: " << mMatch2Count << "(" << double(mMatch2Count)/mTotalToMatch << ")" << std::endl;
		ostr << "Large error count: " << mLargeErrorCount << "(" << double(mLargeErrorCount)/mTotalToMatch << ")" << std::endl;
		ostr << "Average error by edge length: " << mAvgErrInEdgeLength << std::endl;
	}
	if (withRelativeError) ostr << "Average relative error: " << mAvgRelativeError << std::endl;

	ostr << "******************************" << std::endl;
}
