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
const double DiffusionShapeMatcher::DEFAULT_MATCH_THRESH        = 0.52;
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

double distFeature(const HKSFeature& hf1, const HKSFeature& hf2, const VectorND& sig1, const VectorND& sig2, double& tl, int& tn)
{
	tl = max(hf1.m_tl, hf2.m_tl);	// now both are 10.0
	double tu = min(hf1.m_tu, hf2.m_tu);
	tn = int(std::log(tu/tl)/std::log(2.0) + 1.5);	//only consider the overlapping times

	if(tn == 0) return 1.0;

	double dist = 0.0;
	for(int i = 0; i < tn; i++)
	{
		dist += std::pow(sig1.m_vec[i] - sig2.m_vec[i], 2);
	}
	return dist/tn;
}

double distFeaturePair(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const MatchPair& mp1, const MatchPair& mp2)
{
	double t = max(mp1.m_tl, mp2.m_tl); // take the bigger one
	int n = min(mp1.m_tn, mp2.m_tn);

	int x1 = mp1.m_idx1;
	int x2 = mp1.m_idx2;
	int y1 = mp2.m_idx1;
	int y2 = mp2.m_idx2;

	VectorND v1(n), v2(n);
	for(int i = 0; i < n; i++)
	{
		v1.m_vec[i] =pmp1->getVertexPairHK(x1, y1, t);
		v2.m_vec[i] = pmp2->getVertexPairHK(x2, y2, t);
		t *= 2.0;
	}

	return v1.calDistance2(v2) / (v1.length2() + v2.length2());
}

struct AffCand
{
	int idx;
	double val;
	AffCand(int i, double v) : idx(i), val(v) {}
	friend bool operator > (const AffCand& af1, const AffCand& af2)
	{
		if (af1.val > af2.val) return true;
		else return false;
	}
	friend bool operator < (const AffCand& af1, const AffCand& af2)
	{
		if (af1.val > af2.val) return true;
		else return false;
	}
};

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
	if (level == -1)
		return matchedPairsFine;
	else
	{
		//TODO: return match pairs on different levels
		return matchedPairsFine;
	}
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

void DiffusionShapeMatcher::matchFeatures( double matchThresh )
{
	ofstream ostr("output/FeatureMatch.log", ios::trunc);

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
	
////////////////    1. Match coarse features    ////////////////
//////////////////////////////////////////////////////////////////////////
	int size1 = (int) vftCoarse1.size();
	int size2 = (int) vftCoarse2.size();
	vector<VectorND> vsig1(size1), vsig2(size2);
	
	for(int i1 = 0; i1 < size1; i1++)
		calVertexSignature(pOriginalProcessor[0], vftCoarse1[i1], vsig1[i1]);
			
	for(int i2 = 0; i2 < size2; i2++)
		calVertexSignature(pOriginalProcessor[1], vftCoarse2[i2], vsig2[i2]);

 	ostr << "  Coarse features 1: ";
 	for (int i = 0; i < size1; ++i)
 		ostr << vftCoarse1[i].m_index << ' ';
 	ostr << "\n Coarse features 2: ";
 	for (int i = 0; i < size2; ++i)
 		ostr << vftCoarse2[i].m_index << ' ';
 	ostr << endl;
	
	ostr << "candidates: ";
	double sigma1 = 4.0 * matchThresh;
	for(int i = 0; i < size1; i++)
	{
		for(int j = 0; j < size2; j++)
		{
			if (vftCoarse1[i].m_scale != vftCoarse2[j].m_scale) continue;
			double tl = 0.0;
			int tn = 0;
			double d = distFeature(vftCoarse1[i], vftCoarse2[j], vsig1[i], vsig2[j], tl, tn);	//average on each overlapped scale
			if(d < matchThresh)
			{
				ostr << '(' << vftCoarse1[i].m_index <<','<< vftCoarse2[j].m_index << ") "; 

				double score = std::exp(-d/sigma1);
				vTmpMatchPairs.push_back(MatchPair(vftCoarse1[i].m_index, vftCoarse2[j].m_index, tl, tn, score));
				vFeatureMatchScores.push_back(score);
			}
		}
	}
	ostr << endl;

/// create affinity matrix (compatibility of each candidate match)
	const int affinitySize = (int)vTmpMatchPairs.size();
    ostr << "\n  Affinity Matrix size: " << affinitySize << endl;

	mxArray *AM, *VM, *VA;
	AM = mxCreateDoubleMatrix(affinitySize, affinitySize, mxREAL);
	double *am = mxGetPr(AM);
	double sigma2 = 0.02;//0.02; 0.1

	for(int i = 0; i < affinitySize; i++)
	{
		am[i*affinitySize+i] = vFeatureMatchScores[i]; // diagonal

		for(int j = i+1; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i].m_idx1 == vTmpMatchPairs[j].m_idx1)^(vTmpMatchPairs[i].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor conflict; 1-to-multiple not desired? 
			{
				am[i*affinitySize+j] = am[j*affinitySize+i] = 0.0;
				continue;
			}
			double ds = 0.0;
			ds = distFeaturePair(pOriginalProcessor[0], pOriginalProcessor[1], vTmpMatchPairs[i], vTmpMatchPairs[j]);
			am[i*affinitySize+j] = am[j*affinitySize+i] = exp(-ds/sigma2);
		}
	}

//  cout << "affinity matrix done" << endl;
	// solving the greatest eigen-vector
	engPutVariable(m_ep, "AM", AM);
	engEvalString(m_ep, "[VM,VA] = spectral_embedding(AM);");	//computing leading eigenvector of A using svd
	VM = engGetVariable(m_ep, "VM");
	double *vm = mxGetPr(VM);
	VA = engGetVariable(m_ep, "VA");
	double *va = mxGetPr(VA);		//?? not referenced?

	const double c_thresh = 0.01;	//was 0.1
	std::vector<MatchPair> mpc1;

	while(1)
	{
		int i_max = -1;
		double v_max = 0;
		for(int i = 0; i < affinitySize; i++)
		{
			if(abs(vm[i]) > v_max)
			{
				v_max = abs(vm[i]);
				i_max = i;
			}
		}

		if (v_max <= c_thresh) 
		{
			ostr << "  Discarded v_max: " << v_max << endl;
			break;
		}

		ostr << mpc1.size() << ": " << i_max << ',' << v_max << endl;

		bool hit = true;
		const int curMatchSize = (int)mpc1.size();
//*
		if (curMatchSize >= 4)
		{
			int refIndex[4];
			for (int i = 0; i < 4; ++i)
				refIndex[i] =  int( double(curMatchSize)/4.0 * (i + (double)rand()/(double)RAND_MAX)); 
			
			double errorSum(0);
			for (int i = 0; i < 4; ++i)
			{
				int selectIdx = refIndex[i];
				const MatchPair &mp1 = mpc1[selectIdx], &mp2 = vTmpMatchPairs[i_max];
				double geodist1 = mesh1->getGeodesic(mp1.m_idx1, mp2.m_idx1),
					   geodist2 = mesh2->getGeodesic(mp1.m_idx2, mp2.m_idx2);
				double distError = std::abs(geodist1 - geodist2);

				if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 3*mesh1->getAvgEdgeLength() : 0.3*geodist1) )
				{
					hit = false;
					break;
				}
				errorSum += distError;
			}
			if (hit)
			{
				vector<AffCand> mpCandidates;
				for (int j = 0; j < affinitySize; ++j)
				{
					if (j != i_max && 
						vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1 && 
						vTmpMatchPairs[i_max].m_idx2 != vTmpMatchPairs[j].m_idx2 &&
						abs(vm[j]) > c_thresh
						)
					{
						mpCandidates.push_back(AffCand(j, abs(vm[j])));
					}
				}

				int candSize = (int)mpCandidates.size();

				for (int k = 0; k < candSize; ++k)
				{
					double errorSumK(0);
					bool pass = true;
					for (int i = 0; i < 4; ++i)
					{
						int selectIdx = refIndex[i];

						const MatchPair &mp1 = mpc1[selectIdx], &mp2 = vTmpMatchPairs[mpCandidates[k].idx];
						double geodist1 = mesh1->getGeodesic(mp1.m_idx1, mp2.m_idx1),
							   geodist2 = mesh2->getGeodesic(mp1.m_idx2, mp2.m_idx2);
						double distError = std::abs(geodist1 - geodist2);

						if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 3*mesh1->getAvgEdgeLength() : 0.3*geodist1) )
						{
							pass = false;
							break;
						}
						errorSumK += distError;
					}
					if (pass && errorSumK < errorSum)
					{
						i_max = mpCandidates[k].idx;
						errorSum = errorSumK;
					}
				}
			}
		}
//*/
		vm[i_max] = 0.0;
		if (!hit) continue;
		// hit
		mpc1.push_back(vTmpMatchPairs[i_max]);
		
		// now a max candidate is found, remove conflicting candidates
		for(int j = 0; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1) ^ (vTmpMatchPairs[i_max].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor, conflict
				vm[j] = 0.0;
		}
	}// end of while()

	matchedPairsCoarse = mpc1;

	cout << "  Coarse Match computed!" << endl;

	mxDestroyArray(VM);
	mxDestroyArray(VA);
	mxDestroyArray(AM);

	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); /*++iter*/ )
	{
		if (iter->m_idx1 == iter->m_idx2)
		{
			iter++;
			continue;
		}

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
		{
			iter->m_note = -1;
			iter = matchedPairsCoarse.erase(iter);
			//			iter++;
		}
		else iter++;
	}

//////////////////////////////////////////////////////////////////////////
////////    2. Match rest features if possible    ////
	cout << "Mesh 1 features: " << endl;
	for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
	{
		cout << "  "<< hiter->m_index  << ',' << hiter->m_scale;
	}
	cout << endl;
	cout << "Mesh 2 features: " << endl;
	for (vector<HKSFeature>::const_iterator hiter = vftFine2.begin(); hiter != vftFine2.end(); ++hiter)
	{
		cout << "  " << hiter->m_index  << ',' << hiter->m_scale;
	}
	cout << endl;	
	
	matchedPairsFine = matchedPairsCoarse;

	const int coarseMatchSize = (int) matchedPairsCoarse.size();
	vector<int> coarseFeat1, coarseFeat2;
	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); ++iter)
	{
		coarseFeat1.push_back(iter->m_idx1);
		coarseFeat2.push_back(iter->m_idx2);
	}
	vector<int> restFeat1, restFeat2;
	map<int, vector<HKSFeature>::const_iterator> restIdxToHKS1, restIdxToHKS2;
	for (vector<HKSFeature>::const_iterator iter = vftFine1.begin(); iter != vftFine1.end(); ++iter)
	{
		if (find(coarseFeat1.begin(), coarseFeat1.end(), iter->m_index) == coarseFeat1.end())
		{
			restFeat1.push_back(iter->m_index);
			restIdxToHKS1[iter->m_index] = iter;
		}
	}
	for (vector<HKSFeature>::const_iterator iter = vftFine2.begin(); iter != vftFine2.end(); ++iter)
	{
		if (find(coarseFeat2.begin(), coarseFeat2.end(), iter->m_index) == coarseFeat2.end())
		{
			restFeat2.push_back(iter->m_index);
			restIdxToHKS2[iter->m_index] = iter;
		}
	}

	const int restSize1 = (int) restFeat1.size(),
		      restSize2 = (int) restFeat2.size();
	cout << "  Rest Size1: " << restSize1 << '\t' << "Rest Size2: " << restSize2 << endl;
	vector<VectorND> vHKC1, vHKC2;
	vHKC1.resize(restSize1);
	vHKC2.resize(restSize2);
	for (int j = 0; j < restSize1; ++j)
	{
		VectorND& vnd = vHKC1[j];
		vnd.reserve(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd.m_vec[i] = pOriginalProcessor[0]->getVertexPairHK(restFeat1[j], matchedPairsCoarse[i].m_idx1, 80);
//			vnd.m_vec[i] = mesh1->CalGeodesic(restFeat1[j], matchedPairsCoarse[i].m_idx1);
		}
	}
	for (int j = 0; j < restSize2; ++j)
	{
		VectorND& vnd = vHKC2[j];
		vnd.reserve(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd.m_vec[i] = pOriginalProcessor[1]->getVertexPairHK(restFeat2[j], matchedPairsCoarse[i].m_idx2, 80);
//			vnd.m_vec[i] = mesh2->CalGeodesic(restFeat2[j], matchedPairsCoarse[i].m_idx2);
		}
	}

	vector<MatchPair> restMatches;
	for (int i = 0; i < restSize1; ++i)
	{
		for (int j = 0; j < restSize2; ++j)
		{
			if (vftFine2[restFeat1[i]].m_scale != vftFine2[restFeat2[j]].m_scale)
				continue;
			double s = std::exp(-vHKC1[i].calDistance(vHKC2[j]));
			restMatches.push_back(MatchPair(restFeat1[i], restFeat2[j], s));
		}
	}
	cout << "  candidate rest matches: " << restMatches.size() << endl;
	while (!restMatches.empty())
	{
		vector<MatchPair>::iterator maxIter = restMatches.end();
		double maxVal = -1;
		for (vector<MatchPair>::iterator iter = restMatches.begin(); iter != restMatches.end(); ++iter)
		{
			if (iter->m_score > maxVal)
			{
				maxVal = iter->m_score;
				maxIter = iter;
			}
		}
		int maxIdx1 = maxIter->m_idx1, maxIdx2 = maxIter->m_idx2;

		bool pass = true;
		int matchSize = (int) matchedPairsFine.size();
		if (matchSize >= 5)	// from already registered features to judge the validity of new matches
		{
			for (int k = 0; k < 5; ++k)
			{
				int selectIdx = int( double(matchSize)/5.0 * (k + double(rand()/double(RAND_MAX)))); 
				//				int selectIdx = rand() % matchSize;
				const MatchPair &mp1 = matchedPairsFine[selectIdx];
				double geodist1 = mesh1->getGeodesic(mp1.m_idx1, maxIdx1),
					   geodist2 = mesh2->getGeodesic(mp1.m_idx2, maxIdx2);
				double distError = std::abs(geodist1 - geodist2);

				if ( distError >= (geodist1 < 10*mesh1->getAvgEdgeLength() ? 2*mesh1->getAvgEdgeLength() : 0.2*geodist1) ) 
				{
					pass = false;
					break;
				}
			}
		}
		if (pass) 
		{
			matchedPairsFine.push_back(MatchPair(maxIdx1, maxIdx2));
			cout << "  add new match: " << maxIdx1 << ", " << maxIdx2 << endl;
		}
		else 
		{
			restMatches.erase(maxIter);
			continue;
		}

		for (vector<MatchPair>::iterator iter = restMatches.begin(); iter != restMatches.end(); )
		{
			if (iter->m_idx1 == maxIdx1 || iter->m_idx2 == maxIdx2)
				restMatches.erase(iter);
			else iter++;
		}
	}

	ostr << "  Fine matches computed" << endl;
	ostr << " --- Total Matched Features: " << matchedPairsFine.size() << endl;
	ostr << " --- Mesh 1 features num: " << vftFine1.size() << endl;
	ostr << " --- Mesh 2 features num: " << vftFine2.size() << endl;
//// add note to matchings with great discrepancy
	
	for (vector<MatchPair>::iterator iter = matchedPairsFine.begin(); iter != matchedPairsFine.end(); /*++iter*/ )
	{
		for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
		{
			if (hiter->m_scale == 3 && hiter->m_index == iter->m_idx1)
				cout << " --- Scale-3 Feature Pair: " << iter->m_idx1 << ", " << iter->m_idx2 << endl;
		}

		if (iter->m_idx1 == iter->m_idx2)
		{
			iter++;
			continue;
		}

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
		{
			iter->m_note = -1;
			iter = matchedPairsFine.erase(iter);
//			iter++;
		}
		else iter++;
	}
	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); ++iter )
	{
		if (iter->m_idx1 == iter->m_idx2) 
			continue;

		if (!mesh2->isInNeighborRing(iter->m_idx1, iter->m_idx2, 5))	
			iter->m_note = -1;
		//iter = matchedPairs.erase(iter);
	}
}

void DiffusionShapeMatcher::calVertexSignature( const DifferentialMeshProcessor* pOriginalProcessor, const HKSFeature& hf, VectorND& sig) const
{
	double t = hf.m_tl;
	sig.reserve(hf.m_tn);
	for(int i = 0; i < hf.m_tn; i++)
	{
		double eref = 4.0*PI*t;
		sig.m_vec[i] = log(pOriginalProcessor->getVertexHKS(hf.m_index, t) * eref);	//normalization to balance HKS at different t
		t *= 2.0;
	}
}






