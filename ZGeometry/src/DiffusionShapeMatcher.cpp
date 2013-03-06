#include <fstream>
#include <sstream>
#include <exception>
#include <set>
#include <algorithm>
#include <cassert>
#include "DiffusionShapeMatcher.h"
#include "OutputHelper.h"

#define INFINITY 1e10
#define LOCAL_ANCHORS_NUM 8
using namespace std;

extern OutputHelper qout;
extern QString qformat;

const double DiffusionShapeMatcher::DEFAULT_C_RATIO				= 0.2;
const double DiffusionShapeMatcher::DEFAULT_RANK_EPSILON		= 1e-4;
const double DiffusionShapeMatcher::SPARSIFY_EPSILON			= 1e-6;
const double DiffusionShapeMatcher::DEFAULT_FEATURE_TIMESCALE	= 30.0;
const double DiffusionShapeMatcher::DEFAULT_T_MULTIPLIER		= 3.0;
const double DiffusionShapeMatcher::DEFAULT_MATCH_TIME_LOW		= 10.0;
const double DiffusionShapeMatcher::DEFAULT_REGISTER_TIMESCALE	= 80;
const int	 DiffusionShapeMatcher::DEFAULT_NBASE				= 200;//300;
const int	 DiffusionShapeMatcher::NUM_OF_EIGVAL_FOR_ESTIMATE	= 50;
const int	 DiffusionShapeMatcher::DEFAULT_PYRAMID_LEVELS		= 3;
const int	 DiffusionShapeMatcher::MAXIMAL_PYRAMID_LEVELS		= 5;

DiffusionShapeMatcher::DiffusionShapeMatcher()
{
	m_ep = NULL;
	pOriginalProcessor[0] = pOriginalProcessor[1] = NULL; 
	pOriginalMesh[0] = pOriginalMesh[1] = NULL;
	m_bPyramidBuilt = false;
	m_nCurrentMatchLevel = -1;
}

DiffusionShapeMatcher::~DiffusionShapeMatcher()
{
	for (auto iter = liteMP[0].begin(); iter != liteMP[0].end(); ++iter)
	{
		delete *iter;
	}
	for (auto iter = liteMP[1].begin(); iter != liteMP[1].end(); ++iter)
	{
		delete *iter;
	}
}

void DiffusionShapeMatcher::initialize( DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep )
{
	m_ep = ep;
	pOriginalProcessor[0] = pMP1;
	pOriginalProcessor[1] = pMP2;
	pOriginalMesh[0] = pMP1->getMesh();
	pOriginalMesh[1] = pMP2->getMesh();

	liteMP[0].push_back(pMP1);
	liteMP[1].push_back(pMP2);

	meshPyramids[0].setInitialMesh(pOriginalMesh[0]);
	meshPyramids[1].setInitialMesh(pOriginalMesh[1]);
}

CMesh* DiffusionShapeMatcher::getMesh( int obj, int level /*= 0*/ ) const
{
	if (level == 0)
		return this->pOriginalMesh[obj];
	else
		return meshPyramids[obj].getMesh(level);
}

std::vector<MatchPair> DiffusionShapeMatcher::getFeatureMatches( int level ) const
{
	vector<MatchPair> retv;
	return retv;
}

void DiffusionShapeMatcher::constructPyramid( int n )
{
	assert(n >= 1);

	m_nPyramidLevels = n;

	meshPyramids[0].setLevel(n);
	meshPyramids[1].setLevel(n);

	meshPyramids[0].construct();
	meshPyramids[1].construct();

	for (int k = 1; k < n; ++k)
	{
		liteMP[0].push_back(new DifferentialMeshProcessor(meshPyramids[0].getMesh(k)));
		liteMP[1].push_back(new DifferentialMeshProcessor(meshPyramids[1].getMesh(k)));
	}

}





