#include <fstream>
#include <sstream>
#include <exception>
#include <set>
#include <algorithm>
#include "DiffusionShapeMatcher.h"

#define INFINITY 1e10
#define LOCAL_ANCHORS_NUM 8
using namespace std;

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
	pMP[0] = pMP[1] = NULL; 
	pOriginalMesh[0] = pOriginalMesh[1] = NULL;

	m_bPyramidBuilt = false;
	m_nCurrentMatchLevel = -1;
}

void DiffusionShapeMatcher::initialize( DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep )
{
	m_ep = ep;
	pMP[0] = pMP1;
	pMP[1] = pMP2;
	pOriginalMesh[0] = pMP[0]->getMesh();
	pOriginalMesh[1] = pMP[1]->getMesh();
}

CMesh* DiffusionShapeMatcher::getMesh( int obj, int level /*= 0*/ ) const
{
	return this->pOriginalMesh[obj];
//	return meshPyramids[obj].getMesh(level);
}

std::vector<MatchPair> DiffusionShapeMatcher::getFeatureMatches( int level ) const
{
	vector<MatchPair> retv;
	return retv;
}



