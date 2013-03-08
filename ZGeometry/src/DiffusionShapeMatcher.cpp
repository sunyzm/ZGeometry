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

double findTmax( const CMesh* tmesh, int s )
{
	if (tmesh->getBoundaryVertexNum() == 0 ) return -1;
	double geo = tmesh->getGeodesicToBoundary(s) / tmesh->getAvgEdgeLength();	//normalized by average edge length; requires isHole initialized
	return geo*geo/4.0;
}

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
	for (auto iter = liteMP[0].begin()+1; iter != liteMP[0].end(); ++iter)
	{
		delete *iter;
	}
	for (auto iter = liteMP[1].begin()+1; iter != liteMP[1].end(); ++iter)
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

class ExtremaPoint
{
public:
	int index;
	double val;
	ExtremaPoint(int i, double v) : index(i), val(v) {}

	friend bool operator<(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val < ep2.val) return true;
		else return false;
	}

	friend bool operator<=(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val <= ep2.val) return true;
		else return false;
	}


	friend bool operator>(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val > ep2.val) return true;
		else return false;
	}

	friend bool operator>=(const ExtremaPoint& ep1, const ExtremaPoint& ep2)
	{
		if (ep1.val >= ep2.val) return true;
		else return false;
	}
};

void DiffusionShapeMatcher::detectFeatures( int obj, int ring /*= 2*/, int scale /*= 1*/, double baseTvalue /*= DEFAULT_FEATURE_TIMESCALE*/, double talpha /*= DEFAULT_T_MULTIPLIER*/, double thresh /*= 0.04*/ )
{
	assert(obj == 0 || obj == 1);

	const CMesh* fineMesh = pOriginalMesh[obj];
	DifferentialMeshProcessor* pMP = pOriginalProcessor[obj];
	const int fineSize = fineMesh->getVerticesNum();

	vector<HKSFeature>& vF = vFeatures[obj];
	vF.clear();

	for (int s = scale-1; s >= 0; --s)
	{
		double tvalue = baseTvalue * std::pow(talpha, s);
		vector<double> hksv;
		pMP->calKernelSignature(tvalue, HEAT_KERNEL, hksv);

		double sref = 4.0 * PI * tvalue;
		for_each (hksv.begin(), hksv.end(), [=](double& v){ v = std::log(v * sref); });		

		vector<int> vFeatureIdx;
		fineMesh->extractExtrema(hksv, ring, thresh, vFeatureIdx);

		priority_queue<ExtremaPoint, vector<ExtremaPoint>, greater<ExtremaPoint> > extremaPointsQueue;
		for (auto iter = vFeatureIdx.begin(); iter != vFeatureIdx.end(); ++iter)
			extremaPointsQueue.push(ExtremaPoint(*iter, std::abs(hksv[*iter])));
		
		while (!extremaPointsQueue.empty())
		{
			ExtremaPoint ep = extremaPointsQueue.top();
			extremaPointsQueue.pop();

			bool notCoveredFine = true;
			for (auto iter = vF.begin(); iter != vF.end(); ++iter)
			{
				if (fineMesh->isInNeighborRing(iter->m_index, ep.index, 5))
				{
					notCoveredFine = false;
					break;
				}
			}
			if (notCoveredFine)
			{
				vF.push_back(HKSFeature(ep.index, scale));
			}
		}
	}

	double tl = DEFAULT_MATCH_TIME_LOW;	//10.0
	for (auto iter = vF.begin(); iter != vF.end(); ++iter)
	{
		const int vi = iter->m_index;
		double tu = findTmax(fineMesh, vi);	//?? max time to boundary? why geo*geo/4?
		int tn;
		if (tu < 0)
		{
			tn = 8;
			tu = std::pow(2.0, tn-1) * tl;				// no boundary
		}
		else
		{
			if(tu < tl) continue;
			tn = (int)ceil(log(tu/tl)/log(2.0));	// tn is the number of timescales
		}

		iter->setTimes(tl, tu, tn);
	}

	MeshFeatureList *mfl = new MeshFeatureList;

	for (auto iter = vF.begin(); iter != vF.end(); ++iter)
	{
		mfl->addFeature(new HKSFeature(*iter));
	}

	mfl->setIDandName(FEATURE_MULTI_HKS, "Feature_multiple_hks");
	pMP->replaceProperty(mfl);
	pMP->setActiveFeaturesByID(FEATURE_MULTI_HKS);

}

void DiffusionShapeMatcher::matchFeatures()
{
	const CMesh *mesh1 = pOriginalMesh[0], *mesh2 = pOriginalMesh[1];
	const vector<HKSFeature>& vftFine1 = vFeatures[0];
	const vector<HKSFeature>& vftFine2 = vFeatures[1];

	vector<HKSFeature> vftCoarse1, vftCoarse2;

	for_each(vftFine1.begin(), vftFine1.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2) vftCoarse1.push_back(f); 
	});
	for_each(vftFine2.begin(), vftFine2.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2) vftCoarse2.push_back(f);
	});

	vector<MatchPair> vTmpMatchPairs;
	vector<double> vFeatureMatchScores;
	

}






