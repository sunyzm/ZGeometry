#include "ShapeMatcher.h"
#include <fstream>
#include <sstream>
#include <exception>
#include <set>
#include <queue>
#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <ppl.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/util.h>
#include <ZGeom/ZGeom.h>
#include "OutputHelper.h"
#include "global.h"

using namespace std;
using ZGeom::VectorPointwiseProduct;
using ZGeom::PI;
using ZGeom::logic_assert;
using ZGeom::runtime_assert;
using ZGeom::VecNd;
using ZGeom::calGeodesic;
using ZGeom::EigenSystem;
using ZGeom::calHeatTrace;

double calHK(MeshHelper* helper, int v1, int v2, double t) {
    return ZGeom::calHeatKernel(helper->getEigenSystem(CotFormula), v1, v2, t);
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
	bool operator()(const MyNote& Left, const MyNote& Right) const { return (Left.m_score < Right.m_score); }
};

typedef std::priority_queue<MyNote, std::vector<MyNote>, NoteCompare> NoteQueue;


const double ShapeMatcher::DEFAULT_C_RATIO				= 0.2;
const double ShapeMatcher::DEFAULT_RANK_EPSILON			= 1e-4;
const double ShapeMatcher::SPARSIFY_EPSILON				= 1e-6;
const double ShapeMatcher::DEFAULT_FEATURE_TIMESCALE	= 30.0;
const double ShapeMatcher::DEFAULT_T_MULTIPLIER			= 3.0;
const double ShapeMatcher::DEFAULT_MATCH_TIME_LOW		= 10.0;
const double ShapeMatcher::DEFAULT_MATCH_THRESH			= 0.52;
const double ShapeMatcher::DEFAULT_EXTREAMA_THRESH		= 0.04;
const double ShapeMatcher::DEFAULT_REGISTER_TIMESCALE	= 80;
const int	 ShapeMatcher::NUM_OF_EIGVAL_FOR_ESTIMATE	= 50;
const int	 ShapeMatcher::DEFAULT_PYRAMID_LEVELS		= 3;
const int	 ShapeMatcher::MAXIMAL_PYRAMID_LEVELS		= 5;

double findTmax(CMesh* tmesh, int s)
{
 	if (!tmesh->hasBoundary() ) return -1;
 	double geo = ZGeom::calGeodesicToBoundary(*tmesh, s) / tmesh->getAvgEdgeLength();	// normalized by average edge length; requires isHole initialized
 	return geo*geo/4.0;
}

void PointParam::clear()
{
	resize(0);
	m_votes = 0.0;
}

PointParam& PointParam::operator= ( const PointParam& hkp )
{
	resize(hkp.size());
	std::copy_n(mVec, mDim, hkp.c_ptr());
	m_votes = hkp.m_votes;
	return *this;
}

void ParamManager::computeHKParam( const std::vector<int>& anchors, double t /*= 30.0*/ )
{
	const int fineSize = pMP->getMesh()->vertCount();
	const int pn = (int)anchors.size();
	para_dim = pn;
	vCoord.resize(fineSize);

	for (int v = 0; v < fineSize; ++v) {
		PointParam& hkp = vCoord[v];
		hkp.resize(pn);
		for (int i = 0; i < pn; ++i)
			hkp[i] = calHK(pMP, v, anchors[i], t);
 
		hkp.m_votes = hkp.norm2();
	}
}

void ParamManager::para_computeHKC( const std::vector<int>& anchors, double t /*= 30.0*/ )
{
	const int fineSize = pMP->getMesh()->vertCount();
	const int pn = (int)anchors.size();
	para_dim = pn;
	vCoord.resize(fineSize);

	Concurrency::parallel_for(0, fineSize, [&](int v){
		PointParam& hkp = vCoord[v];
		hkp.resize(pn);
		for (int i = 0; i < pn; ++i)
			hkp[i] = calHK(pMP, v, anchors[i], t);
		hkp.m_votes = hkp.norm2();
	});
}

void ParamManager::para_computeHKS( const std::vector<double>& times )
{
	const int fineSize = pMP->getMesh()->vertCount();
	const int sn = (int)times.size();
	para_dim = sn;;
	vSignature.resize(fineSize);

	Concurrency::parallel_for(0, fineSize, [&](int v){
		vSignature[v].resize(sn);
		for (int i = 0; i < sn; ++i)
			vSignature[v][i] = calHK(pMP, v, v, times[sn]);
	});
}

double distHKS2(MeshHelper* pmp1, MeshHelper* pmp2, int i1, int i2, double tl, int tn)
{
// 	if (tu < tl) return 1.0;
// 	int tn = int(std::log(tu/tl)/std::log(2.0)) + 1;	//only consider the overlapping times
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

    VecNd v1(tn), v2(tn);
	double t = tl;
	for(int i = 0; i < tn; i++) {
        v1[i] = calHK(pmp1, i1, i1, t) / ZGeom::calHeatTrace(es1, t);
        v2[i] = calHK(pmp2, i2, i2, t) / ZGeom::calHeatTrace(es2, t);
		t *= 2.0;
	}

	return std::pow((v1 - v2).norm2(), 2) / tn;
}

double distHKPair2(MeshHelper* pmp1, MeshHelper* pmp2, const MatchPair& mp1, const MatchPair& mp2, double tl, double tu)
{
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

	if (tu < tl) return 1.0;
	int tn = int(std::log(tu/tl)/std::log(2.0)) + 1;	

	int x1 = mp1.m_idx1;
	int x2 = mp1.m_idx2;
	int y1 = mp2.m_idx1;
	int y2 = mp2.m_idx2;

	VecNd v1(tn), v2(tn);
	double t = tl;
	for(int i = 0; i < tn; i++) {
        v1[i] = calHK(pmp1, x1, y1, t) / ZGeom::calHeatTrace(es1, t);
        v2[i] = calHK(pmp2, x2, y2, t) / ZGeom::calHeatTrace(es2, t);
		t *= 2.0;
	}

	double dist = (v1 - v2).norm2();
	return dist*dist / tn;
}

double distHksFeature2(MeshHelper* pmp1, MeshHelper* pmp2, const HKSFeature& hf1, const HKSFeature& hf2, double& tl, int& tn)
{
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

	tl = max(hf1.m_tl, hf2.m_tl);	// now both are 10.0
	double tu = min(hf1.m_tu, hf2.m_tu);
	if (tu < tl) return 1.0;
	tn = int(std::log(tu/tl)/std::log(2.0)) + 1;	//only consider the overlapping times

	VecNd v1(tn), v2(tn);
	double t = tl;
	for(int i = 0; i < tn; i++) {
        v1[i] = calHK(pmp1, hf1.m_index, hf1.m_index, t) / ZGeom::calHeatTrace(es1, t);
        v2[i] = calHK(pmp2, hf2.m_index, hf2.m_index, t) / ZGeom::calHeatTrace(es2, t);
		t *= 2.0;
	}

	double dist = (v1 - v2).norm2();
	return dist*dist / tn;
}

double distHksFeaturePair2(MeshHelper* pmp1, MeshHelper* pmp2, const MatchPair& mp1, const MatchPair& mp2)
{
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

	double tl = max(mp1.m_tl, mp2.m_tl);	// now both are 10.0
	double tu = min(mp1.m_tu, mp2.m_tu);
	if (tu < tl) return 1.0;
	int tn = int(std::log(tu/tl)/std::log(2.0)) + 1;	//only consider the overlapping times

	int x1 = mp1.m_idx1;
	int x2 = mp1.m_idx2;
	int y1 = mp2.m_idx1;
	int y2 = mp2.m_idx2;

	VecNd v1(tn), v2(tn);
	double t = tl;
	for(int i = 0; i < tn; i++) {
        v1[i] = calHK(pmp1, x1, y1, t) / ZGeom::calHeatTrace(es1, t);
        v2[i] = calHK(pmp2, x2, y2, t) / ZGeom::calHeatTrace(es2, t);
		t *= 2.0;
	}

	double dist = (v1 - v2).norm2();
	return dist * dist / tn;
//	return v1.calDistance2(v2) / (v1.length2() + v2.length2());
}

ShapeMatcher::ShapeMatcher()
{
	m_ep = NULL;
	pOriginalProcessor[0] = pOriginalProcessor[1] = NULL; 
	pOriginalMesh[0] = pOriginalMesh[1] = NULL;
	m_bInitialized = false;
	m_bPyramidBuilt = false;
	m_bHasGroundTruth = false;
	m_nAlreadyRegisteredLevel = -1;

	m_vFeatures.resize(2);
}

ShapeMatcher::~ShapeMatcher()
{
	if (!m_bInitialized) return;

	int mpSize1 = liteMP[0].size(), mpSize2 = liteMP[1].size();
	for (int k = 1; k < mpSize1; ++k) delete liteMP[0][k];
	for (int k = 1; k < mpSize2; ++k) delete liteMP[1][k];
	
	std::cout << "DiffusionShapeMatcher destroyed!" << std::endl;
}

void ShapeMatcher::initialize( MeshHelper* pMP1, MeshHelper* pMP2, Engine *ep )
{
	m_ep = ep;
	pOriginalProcessor[0] = pMP1;
	pOriginalProcessor[1] = pMP2;
	pOriginalMesh[0] = pMP1->getMesh();
	pOriginalMesh[1] = pMP2->getMesh();

	liteMP[0].push_back(pMP1);
	liteMP[1].push_back(pMP2);

	m_ParamMgr[0].initialize(pMP1);
	m_ParamMgr[1].initialize(pMP2);
	
	m_nPyramidLevels = 1;
	setRegistrationLevels(1);
	m_bInitialized = true;
}

CMesh* ShapeMatcher::getMesh( int obj, int level /*= 0*/ ) const
{
	if (level == 0) return this->pOriginalMesh[obj];
    else return nullptr;
//	else return meshPyramids[obj].getMesh(level);
}

void ShapeMatcher::detectFeatures( int obj, int ring /*= 2*/, int nScales /*= 1*/, double baseTvalue /*= DEFAULT_FEATURE_TIMESCALE*/, double talpha /*= DEFAULT_T_MULTIPLIER*/, double thresh /*= DEFAULT_EXTREAMA_THRESH*/ )
{
	logic_assert(obj == 0 || obj == 1);
	m_vFeatures[obj].clear();

	CMesh* fineMesh = pOriginalMesh[obj];
	MeshHelper* pMP = pOriginalProcessor[obj];
	const int fineSize = fineMesh->vertCount();
	vector<HKSFeature>& vF = m_vFeatures[obj];
	vF.clear();

	vector<double> vScaleValues(nScales);
	for (int s = nScales-1; s >= 0; --s) 
		vScaleValues[s] = baseTvalue * std::pow(talpha, s);

	for (int s = nScales-1; s >= 0; --s) {
		
		cout << "=== Detection timescale: " << vScaleValues[s] << " ===" << endl;
        vector<double> hksv = ZGeom::calHeatKernelSignature(pMP->getEigenSystem(SymCot), vScaleValues[s]);
		double sref = 4.0 * PI * vScaleValues[s];
		transform(hksv.begin(), hksv.end(), hksv.begin(), [=](double v){ return std::log(v * sref);} );		
                
        // <index, minOrMax>
        vector<pair<int, int>> vFeatureIdx = ZGeom::extractSignedMeshExtrema(*fineMesh, hksv, ring, thresh);

		for (auto iter = vFeatureIdx.begin(); iter != vFeatureIdx.end(); ++iter) {
			if ( vF.end() == find_if(vF.begin(), vF.end(), 
				//[&](const HKSFeature& feat){return fineMesh->isInNeighborRing(feat.m_index, iter->first, 1);})
				[&](const HKSFeature& feat){return feat.m_index == iter->first;}) )	{
					vF.push_back(HKSFeature(iter->first, s, iter->second)); 
			}			
		}

// 		priority_queue<ExtremaPoint, vector<ExtremaPoint>, greater<ExtremaPoint> > extremaPointsQueue;
// 		for_each (vFeatureIdx.begin(), vFeatureIdx.end(), [&](int idx){
// 			extremaPointsQueue.push(ExtremaPoint(idx, std::abs(hksv[idx])));
// 		});
// 		while (!extremaPointsQueue.empty())
// 		{
// 			ExtremaPoint ep = extremaPointsQueue.top();
// 			extremaPointsQueue.pop();
// 
// 			if ( vF.end() == find_if(vF.begin(), vF.end(), 
// 				[&](const HKSFeature& feat){return fineMesh->isInNeighborRing(feat.m_index, ep.index, 3);})
// 			   )
// 				vF.push_back(HKSFeature(ep.index, s));
// 		}
	}

	/* Begin LOG-FEATURES */
	cout << "[Candidate features " << obj << "]\n";
	std::sort(vF.begin(), vF.end(), [](HKSFeature& f1, HKSFeature& f2) {return f1.m_index < f2.m_index;});
	std::sort(vF.begin(), vF.end(), [](HKSFeature& f1, HKSFeature& f2) {return f1.m_scale > f2.m_scale;});
	for (auto iter = vF.begin(); iter != vF.end(); ++iter)
		cout << "(" << iter->m_index << "," << iter->m_scale << ")" << ", ";
	cout << endl;
	/* End LOG-FEATURES */

	double tl = g_configMgr.getConfigValueDouble("MATCH_TIME_LOWEST");	// 10.0
	int default_num_scales = g_configMgr.getConfigValueInt("NUM_MATCH_SCALES");	// 8
	for (auto iter = vF.begin(); iter != vF.end(); )
	{
		const int vi = iter->m_index;
		double tu = findTmax(fineMesh, vi);	//?? max time to boundary? why geo*geo/4?
		int tn;
		if (tu < 0)	// no boundary
		{
			tn = default_num_scales;
			tu = tl * std::pow(2.0, tn-1);
		} else {
			if(tu < tl) { tu = tl;/*iter = vF.erase(iter); continue;*/}	// features too close to boundary are discarded
			tn = log(tu/tl)/log(2.0) + 1.01;	    // tn is the number of timescales
		}
		iter->setTimes(tl, tu, tn);
		++iter;
	}

	MeshFeatureList *mfl = new MeshFeatureList;
	for (auto iter = vF.begin(); iter != vF.end(); ++iter) mfl->addFeature(new HKSFeature(*iter));
	fineMesh->addAttrMeshFeatures(*mfl, StrAttrFeatureMultiHKS);

	m_bFeatureDetected = true;
}

void ShapeMatcher::matchFeatures( std::ostream& flog, double matchThresh /*= DEFAULT_MATCH_THRESH*/ )
{
	CMesh *mesh1 = pOriginalMesh[0], *mesh2 = pOriginalMesh[1];
	const vector<HKSFeature>& vftFine1 = m_vFeatures[0];
	const vector<HKSFeature>& vftFine2 = m_vFeatures[1];
	
	vector<HKSFeature> vftCoarse1, vftCoarse2;
	std::vector<MatchPair> matchedPairsCoarse, matchedPairsFine;

	flog << "Before defining vftCoarse" << endl;

	for_each(vftFine1.begin(), vftFine1.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2) vftCoarse1.push_back(f); 
	});
	for_each(vftFine2.begin(), vftFine2.end(), [&](const HKSFeature& f){ 
		if(f.m_scale >= 2) vftCoarse2.push_back(f);
	});

	std::sort(vftCoarse1.begin(), vftCoarse1.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });
	std::sort(vftCoarse2.begin(), vftCoarse2.end(), [](const HKSFeature& f1, const HKSFeature& f2) { return f1.m_index < f2.m_index; });

	vector<MatchPair> vTmpMatchPairs;
	vector<double> vFeatureMatchScores;
	
////////////////    1. Match coarse features    ////////////////
//////////////////////////////////////////////////////////////////////////
	int size1 = (int) vftCoarse1.size();
	int size2 = (int) vftCoarse2.size();
// 	vector<VectorND> vsig1(size1), vsig2(size2);
	
	flog << "-- vftCoarse1: " << size1 << "; vftCoarse2: " << size2 << "\nBefore calVertexSignature" << endl;
	flog << "-- Coarse features 1: ";
	for (int i = 0; i < size1; ++i)
		flog << vftCoarse1[i].m_index << ' ';
	flog << "\n-- Coarse features 2: ";
	for (int i = 0; i < size2; ++i)
		flog << vftCoarse2[i].m_index << ' ';
	flog << endl;

// 	vector<HKSFeature>* vftCoarse[2] = {&vftCoarse1, &vftCoarse2};
// 	vector<VectorND>* vsig[2] = {&vsig1, &vsig2};
// 
// 	Concurrency::parallel_for(0, 1, [](int obj) {
// 		int size = vftCoarse[obj]->size();
// 		for(int i = 0; i < size; i++)
// 			calVertexSignature(pOriginalProcessor[obj], vftCoarse[obj]->at(i), vsig[obj]->at(i));
// 	});

#if (0)//(_MSC_VER >= 1600)
	Concurrency::parallel_for (0, size1, [&](int i1) {
		calVertexSignature(pOriginalProcessor[0], vftCoarse1[i1], vsig1[i1]);
	});
	Concurrency::parallel_for(0, size2, [&](int i2) {
		calVertexSignature(pOriginalProcessor[1], vftCoarse2[i2], vsig2[i2]);
	});
#else
//   	for(int i1 = 0; i1 < size1; i1++)
//   		calVertexSignature(pOriginalProcessor[0], vftCoarse1[i1], vsig1[i1]);
//     for(int i2 = 0; i2 < size2; i2++)
//    		calVertexSignature(pOriginalProcessor[1], vftCoarse2[i2], vsig2[i2]);
#endif

	flog << "-- candidates: ";
	double sigma1 = 4.0 * matchThresh;
	for(int i = 0; i < size1; i++)
	{
		for(int j = 0; j < size2; j++)
		{
			if (vftCoarse1[i].m_scale != vftCoarse2[j].m_scale) continue;
			double tl = 0.0;
			int tn = 0;
			double d = distHksFeature2(pOriginalProcessor[0], pOriginalProcessor[1], vftCoarse1[i], vftCoarse2[j], tl, tn);	//average on each overlapped scale
			if(d < matchThresh)
			{
				flog << '(' << vftCoarse1[i].m_index <<','<< vftCoarse2[j].m_index << ") "; 

				double score = std::exp(-d/sigma1);
				vTmpMatchPairs.push_back(MatchPair(vftCoarse1[i].m_index, vftCoarse2[j].m_index, tl, tn, score));
				vFeatureMatchScores.push_back(score);
			}
		}
	}
	flog << endl;

	/* ---- create affinity matrix (compatibility of each candidate match) ---- */
	const int affinitySize = (int)vTmpMatchPairs.size();
	flog << "\nAffinity Matrix size: " << affinitySize << endl;

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
			ds = distHksFeaturePair2(pOriginalProcessor[0], pOriginalProcessor[1], vTmpMatchPairs[i], vTmpMatchPairs[j]);
			am[i*affinitySize+j] = am[j*affinitySize+i] = exp(-ds/sigma2);
		}
	}
	flog << "Affinity matrix constructed!" << endl;

	/* ---- solving the greatest eigenvector (PCA?) ---- */
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
			flog << "-- Discarded v_max: " << v_max << endl;
			break;
		}

		flog << "  " << mpc1.size() << ": " << i_max << ',' << v_max << endl;

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
				double geodist1 = calGeodesic(*mesh1, mp1.m_idx1, mp2.m_idx1),
                       geodist2 = calGeodesic(*mesh2, mp1.m_idx2, mp2.m_idx2);
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
				vector<ExtremaPoint> mpCandidates;
				for (int j = 0; j < affinitySize; ++j)
				{
					if (j != i_max && 
						vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1 && 
						vTmpMatchPairs[i_max].m_idx2 != vTmpMatchPairs[j].m_idx2 &&
						abs(vm[j]) > c_thresh
						)
					{
						mpCandidates.push_back(ExtremaPoint(j, abs(vm[j])));
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

						const MatchPair &mp1 = mpc1[selectIdx], &mp2 = vTmpMatchPairs[mpCandidates[k].index];
						double geodist1 = calGeodesic(*mesh1, mp1.m_idx1, mp2.m_idx1),
							   geodist2 = calGeodesic(*mesh2, mp1.m_idx2, mp2.m_idx2);
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
						i_max = mpCandidates[k].index;
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
	} // end of while()

	matchedPairsCoarse = mpc1;

	flog << "Coarse Match computed!" << endl;

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
	flog << "------------------------------------------" << endl;
	flog << "\nMatch rest features if possible" << endl;
	flog << "  Mesh 1 features: ";
	for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
		flog << hiter->m_index  << ',' << hiter->m_scale << "  ";
	flog << endl;
	flog << "  Mesh 2 features: ";
	for (vector<HKSFeature>::const_iterator hiter = vftFine2.begin(); hiter != vftFine2.end(); ++hiter)
		flog << hiter->m_index  << ',' << hiter->m_scale << "  ";
	flog << endl;	
	
	matchedPairsFine = matchedPairsCoarse;

	const int coarseMatchSize = (int) matchedPairsCoarse.size();
	vector<int> coarseFeat1, coarseFeat2;
	for (vector<MatchPair>::iterator iter = matchedPairsCoarse.begin(); iter != matchedPairsCoarse.end(); ++iter)
	{
		coarseFeat1.push_back(iter->m_idx1);
		coarseFeat2.push_back(iter->m_idx2);
	}

	flog << "!!test ok 1" << endl;

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
	flog << "!!test ok 2" << endl;

	for (vector<HKSFeature>::const_iterator iter = vftFine2.begin(); iter != vftFine2.end(); ++iter)
	{
		if (find(coarseFeat2.begin(), coarseFeat2.end(), iter->m_index) == coarseFeat2.end())
		{
			restFeat2.push_back(iter->m_index);
			restIdxToHKS2[iter->m_index] = iter;
		}
	}

	const int restSize1 = (int) restFeat1.size(), restSize2 = (int) restFeat2.size();
	flog << "Rest Size1: " << restSize1 << '\t' << "Rest Size2: " << restSize2 << endl;

	vector<VecNd> vHKC1, vHKC2;
	vHKC1.resize(restSize1);
	vHKC2.resize(restSize2);
	for (int j = 0; j < restSize1; ++j)
	{
		VecNd& vnd = vHKC1[j];
		vnd.resize(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd[i] = calHK(pOriginalProcessor[0], restFeat1[j], matchedPairsCoarse[i].m_idx1, 80);
		}
	}

	for (int j = 0; j < restSize2; ++j)
	{
		VecNd& vnd = vHKC2[j];
		vnd.resize(coarseMatchSize);
		for (int i = 0; i < coarseMatchSize; ++i)
		{
			vnd[i] = calHK(pOriginalProcessor[1], restFeat2[j], matchedPairsCoarse[i].m_idx2, 80);
		}
	}

	vector<MatchPair> restMatches;
	for (int i = 0; i < restSize1; ++i)
	{
		for (int j = 0; j < restSize2; ++j)
		{
			if (restIdxToHKS1[restFeat1[i]]->m_scale != restIdxToHKS2[restFeat2[j]]->m_scale)
				continue;
			double s = std::exp(-(vHKC1[i]-vHKC2[j]).norm2());
			restMatches.push_back(MatchPair(restFeat1[i], restFeat2[j], s));
		}
	}

	flog << "Candidate rest matches: " << restMatches.size() << endl;
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
				double geodist1 = calGeodesic(*mesh1, mp1.m_idx1, maxIdx1),
					   geodist2 = calGeodesic(*mesh2, mp1.m_idx2, maxIdx2);
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
			flog << "  add new match: " << maxIdx1 << ", " << maxIdx2 << endl;
		}
		else 
		{
			restMatches.erase(maxIter);
			continue;
		}

		for (vector<MatchPair>::const_iterator iter = restMatches.begin(); iter != restMatches.end(); )
		{
			if (iter->m_idx1 == maxIdx1 || iter->m_idx2 == maxIdx2)
				iter = restMatches.erase(iter);
			else iter++;
		}
	}

	flog << "  Fine matches computed" << endl;
	flog << " --- Total Matched Features: " << matchedPairsFine.size() << endl;
	flog << " --- Mesh 1 features num: " << vftFine1.size() << endl;
	flog << " --- Mesh 2 features num: " << vftFine2.size() << endl;

	flog << "\nResults: " << endl;
	for(auto iter = matchedPairsFine.begin(); iter != matchedPairsFine.end(); ++iter)
		flog << "(" << iter->m_idx1 << ", " << iter->m_idx2 << ")" << endl;

	//* add note to matchings with great discrepancy *//	
	for (vector<MatchPair>::iterator iter = matchedPairsFine.begin(); iter != matchedPairsFine.end(); /*++iter*/ )
	{
		for (vector<HKSFeature>::const_iterator hiter = vftFine1.begin(); hiter != vftFine1.end(); ++hiter)
		{
			if (hiter->m_scale == 3 && hiter->m_index == iter->m_idx1)
				flog << " --- Scale-3 Feature Pair: " << iter->m_idx1 << ", " << iter->m_idx2 << endl;
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

	forceInitialAnchors(matchedPairsFine);
} // DiffusionShapmeMatcher::matchFeatures()

void ShapeMatcher::calVertexSignature(MeshHelper* pOriginalProcessor, const HKSFeature& hf, ZGeom::VecNd& sig)
{
	double t = hf.m_tl;
	sig.resize(hf.m_tn);
	for(int i = 0; i < hf.m_tn; i++)
	{
		double eref = 4.0*PI*t;
		sig[i] = std::log(calHK(pOriginalProcessor, hf.m_index, hf.m_index, t) * eref);	//normalization to balance HKS at different t
		t *= 2.0;
	}
}

void ShapeMatcher::setRegistrationLevels( int val )
{
	m_nRegistrationLevels = min(val, m_nPyramidLevels);
	m_vFeatureMatchingResults.resize(m_nRegistrationLevels + 1);
	m_vRegistrationResutls.resize(m_nRegistrationLevels);

	m_nAlreadyRegisteredLevel = m_nRegistrationLevels - 1;
	m_nAlreadyMatchedLevel = m_nRegistrationLevels;
}

const std::vector<MatchPair>& ShapeMatcher::getMatchedFeaturesResults ( int level ) const
{
	if (level < -1 || level > m_nRegistrationLevels) 
		throw logic_error("Selected level out of bound!");
	else if (level == -1) return m_vFeatureMatchingResults.back();
	else return m_vFeatureMatchingResults[level];	
}

const std::vector<MatchPair>& ShapeMatcher::getRegistrationResults( int level ) const
{
	if (level < -1 || level >= m_nRegistrationLevels) 
		throw logic_error("Selected level out of bound!");
	else if (level == -1) return m_vRegistrationResutls.back();
	else return m_vRegistrationResutls[level];	
}

double ShapeMatcher::computeMatchScore( int idx1, int idx2, double sigma /*= 0.02*/ ) const
{
	const PointParam &hkp1 = m_ParamMgr[0].vCoord[idx1],
					 &hkp2 = m_ParamMgr[1].vCoord[idx2];
	double dist = (hkp1 - hkp2).norm2();
	double dist2 = dist*dist / hkp1.size();
	return std::exp(-dist2 / sigma);
}

void ShapeMatcher::ComputeTensorFeature3(MeshHelper* pmp, int i, int j, int k, double t, double* sang)
{

	if(i==j || i==k || j==k)		//not a triangle
	{
		sang[0] = -10.0;
		sang[1] = -10.0;
		sang[2] = -10.0;
		return;
	}

	double d1 = calHK(pmp, i, j, t);
	double d2 = calHK(pmp, j, k, t);
	double d3 = calHK(pmp, k, i, t);
	if(d1<0.0) d1 = 1e-15;
	if(d2<0.0) d2 = 1e-15;
	if(d3<0.0) d3 = 1e-15;
	d1 = sqrt(-4.0*t*log(d1));
	d2 = sqrt(-4.0*t*log(d2));
	d3 = sqrt(-4.0*t*log(d3));

	// cosine

	sang[0] = (d1*d1+d2*d2-d3*d3)/(2.0*d1*d2);
	sang[1] = (d2*d2+d3*d3-d1*d1)/(2.0*d2*d3);
	sang[2] = (d3*d3+d1*d1-d2*d2)/(2.0*d3*d1);

// sine theorem: area = 2*a*b*sin(C)
//	double s = (d1+d2+d3)/2.0;
//	double a = sqrt(s*(s-d1)*(s-d2)*(s-d3));

//	sang[0] = a/(2*d1*d2);
//	sang[1] = a/(2*d2*d3);
//	sang[2] = a/(2*d3*d1);
}

void ShapeMatcher::ComputeTensorFeature6(MeshHelper* pmp, int i, int j, int k, double t, double* sang, bool sweep /*= false*/ )
{
//  	if(i==j || i==k || j==k)		//not a triangle
//  	{
//  		sang[0] = -100.0;
//  		sang[1] = -100.0;
//  		sang[2] = -100.0;
// 		    sang[3] = sang[4] = sang[5] = -100.;
//  		return;
//  	}
	
    double d1 = calHK(pmp, i, j, t);
    double d2 = calHK(pmp, j, k, t);
    double d3 = calHK(pmp, k, i, t);

    double s1 = calHK(pmp, i, i, t);
    double s2 = calHK(pmp, j, j, t);
    double s3 = calHK(pmp, k, k, t);
	
	sang[0] = d1;
	sang[1] = d2;
	sang[2] = d3;
	sang[3] = s1;
	sang[4] = s2;
	sang[5] = s3;
}


void ShapeMatcher::PairGraphMatching(Engine *ep, MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<HKSFeature>& vFeatures1, const std::vector<HKSFeature>& vFeatures2, std::vector<MatchPair>& vMatchedPair, double para_thresh, bool verbose /*= false*/)
{
	int size1 = (int) vFeatures1.size();
	int size2 = (int) vFeatures2.size();

	vector<MatchPair> vTmpMatchPairs;

	double sigma1 = 4.0 * para_thresh;
	double sigma2 = 0.02;	//0.02; 0.1

	for(int i = 0; i < size1; i++)
	{
		for(int j = 0; j < size2; j++)
		{
			double tl = 0.0;
			int tn = 0;
			double d = distHksFeature2(pmp1, pmp2, vFeatures1[i], vFeatures2[j], tl, tn);	//average on each overlapped scale
			if(d < para_thresh)
			{
				double score = std::exp(-d/sigma1);
				vTmpMatchPairs.push_back(MatchPair(vFeatures1[i].m_index, vFeatures2[j].m_index, tl, tn, score));
			}
		}
	}

	/* ---- create affinity matrix (compatibility of each candidate match) ---- */
	const int affinitySize = (int)vTmpMatchPairs.size();
	mxArray *AM, *VM, *VA;
	AM = mxCreateDoubleMatrix(affinitySize, affinitySize, mxREAL);
	double *am = mxGetPr(AM);
	for(int i = 0; i < affinitySize; i++)
	{
		am[i*affinitySize+i] = vTmpMatchPairs[i].m_score; // diagonal

		for(int j = i+1; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i].m_idx1 == vTmpMatchPairs[j].m_idx1)^(vTmpMatchPairs[i].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor conflict; 1-to-multiple not desired? 
			{
				am[i*affinitySize+j] = am[j*affinitySize+i] = 0.0;
				continue;
			}
			double ds = distHksFeaturePair2(pmp1, pmp2, vTmpMatchPairs[i], vTmpMatchPairs[j]);
			am[i*affinitySize+j] = am[j*affinitySize+i] = exp(-ds/sigma2);
		}
	}

	/* ---- solving the greatest eigenvector (PCA?) ---- */
	engPutVariable(ep, "AM", AM);
	engEvalString(ep, "[VM,VA] = spectral_embedding(AM);");	//computing leading eigenvector of A using svd
	VM = engGetVariable(ep, "VM");
	double *vm = mxGetPr(VM);
	VA = engGetVariable(ep, "VA");
	double *va = mxGetPr(VA);		//?? not referenced?

	const double c_thresh = 0.01;	//was 0.1
	std::vector<MatchPair> mpc1;

	transform(vm, vm+affinitySize, vm, [](double v){ return abs(v); });

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
			if (verbose) cout << "-- Discarded v_max: " << v_max << endl;
			break;
		}

		const int curMatchSize = (int)mpc1.size();		
		if (verbose) cout << "\t" << curMatchSize << ": " << i_max << ',' << v_max << endl;
		
		vm[i_max] = 0.0;
		mpc1.push_back(vTmpMatchPairs[i_max]);
		
		// now a max candidate is found, remove conflicting candidates
		for(int j = 0; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1) ^ (vTmpMatchPairs[i_max].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor, conflict
				vm[j] = 0.0;
		}
	} // end of while()

	mxDestroyArray(VM);
	mxDestroyArray(VA);
	mxDestroyArray(AM);

	vMatchedPair = mpc1;
	if (verbose) cout << "Pair graph matching computed!" << endl;
}

double ShapeMatcher::TensorGraphMatching6( Engine *ep, 
	    MeshHelper* pmp1, MeshHelper* pmp2, 
	    const std::vector<int>& vFeatures1, 
        const std::vector<int>& vFeatures2, 
        std::vector<MatchPair>& matched, double para_t, 
        double para_thresh, bool verbose/* = false */)
{
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

	matched.clear();
	const int vsize1 = (int)vFeatures1.size();	// input feature size 1
	const int vsize2 = (int)vFeatures2.size();    // input feature size 2
	if (vsize1 < 3 || vsize2 < 3) return 0;

	if (vsize1 > vsize2)
	{
		vector<MatchPair> oppositeMatched;
		double score = TensorGraphMatching6(ep, pmp2, pmp1, vFeatures2, vFeatures1, oppositeMatched, para_t, para_thresh, verbose);
		matched.resize(oppositeMatched.size());
		transform(oppositeMatched.begin(), oppositeMatched.end(), matched.begin(), [](const MatchPair& mp){return MatchPair(mp.m_idx2, mp.m_idx1, mp.m_score);});
		return score;
	}

	// generate triangles
	vector<int> triangs;
	for (int i = 0; i < vsize1; i++)
	{
		for (int j = i+1; j < vsize1; j++)
		{
			double c1 = calHK(pmp1, vFeatures1[i], vFeatures1[j], para_t);
			for (int k = j+1; k < vsize1; k++)
			{
				double c2 = calHK(pmp1, vFeatures1[j], vFeatures1[k], para_t),
					   c3 = calHK(pmp1, vFeatures1[i], vFeatures1[k], para_t);
				int count_valid = 0;
				if (c1 >= 1e-5) count_valid++;
				if (c2 >= 1e-5) count_valid++;
				if (c3 >= 1e-5) count_valid++;

				if (count_valid == 0) continue;

				triangs.push_back(i);
				triangs.push_back(j);
				triangs.push_back(k);
			}			
		}
	}
	
	if (triangs.size() < 3)
	{
		cout << "no valid triangles!" << endl;
		return 0;
	}

	int tsize1 = (int)triangs.size() / 3;

	if (verbose)
		cout << "#Total query triangles: " << tsize1 << endl;
	int tsize2 = vsize2 * vsize2 * vsize2;

	// compute feature descriptors
	const int FeatureDim = 6;
	mxArray *mxfeat1, *mxfeat2, *mxtris, *mxnumbs, *mX2, *vX2, *score;
	mxfeat1 = mxCreateDoubleMatrix(FeatureDim, tsize1, mxREAL);
	double *pfeat1 = mxGetPr(mxfeat1);
	mxfeat2 = mxCreateDoubleMatrix(FeatureDim, tsize2, mxREAL);
	double *pfeat2 = mxGetPr(mxfeat2);
	mxtris = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *ptris = mxGetPr(mxtris);
	mxnumbs = mxCreateDoubleMatrix(1, 4, mxREAL);
	double *pnumbs = mxGetPr(mxnumbs);

	const int maxNN = 32;
	pnumbs[0] = vsize1; // nP1
	pnumbs[1] = vsize2; // nP2
	pnumbs[2] = tsize1; // nT
	if(tsize1 > 2*maxNN) pnumbs[3] = maxNN;     // nNN
	else pnumbs[3] = max(1, tsize1 / 2);

	copy(triangs.begin(), triangs.end(), ptris);


	//	for(int i = 0; i < tsize1; i++)
	Concurrency::parallel_for(0, tsize1, [&](int i)
	{
		int vi = triangs[i*3];
		int vj = triangs[i*3+1];
		int vk = triangs[i*3+2];
		ComputeTensorFeature6(pmp1, vFeatures1[vi], vFeatures1[vj], vFeatures1[vk], para_t, pfeat1 + i*FeatureDim);
	}
	);
	Concurrency::parallel_for(0, vsize2, [&](int i)
		//	for(int i = 0; i < vsize2; i++)
	{
		for(int j = 0; j < vsize2; j++)
		{
			for(int k = 0; k < vsize2; k++)
			{
				ComputeTensorFeature6(pmp2, vFeatures2[i], vFeatures2[j], vFeatures2[k], para_t, pfeat2 + ((i*vsize2+j)*vsize2+k)*FeatureDim, /*sweep=*/true);
			}
		}
	}
	);
	// invoke matlab for tensor matching

    double ht1 = calHeatTrace(es1, para_t), ht2 = calHeatTrace(es2, para_t);
	transform(pfeat1, pfeat1+FeatureDim*tsize1, pfeat1, [=](double v){ return v/ht1; });
	transform(pfeat2, pfeat2+FeatureDim*tsize2, pfeat2, [=](double v){ return v/ht2; });

	engPutVariable(ep, "feat1", mxfeat1);
	engPutVariable(ep, "feat2", mxfeat2);
	engPutVariable(ep, "tris", mxtris);
	engPutVariable(ep, "numbs", mxnumbs);
	//	system("PAUSE");	// for manually testing tensorMat
	engEvalString(ep, "[vX2,mX2,score]=tensorMat(feat1,feat2,tris,numbs);");
	mX2 = engGetVariable(ep, "mX2");
	double *px = mxGetPr(mX2);  // best match of each shape-1 feature in the set of shape-2 features
	vX2 = engGetVariable(ep, "vX2");
	double *pv = mxGetPr(vX2);	// corresponding match scores
	score = engGetVariable(ep, "score");	// general matching score
	double *ps = mxGetPr(score);

	/* interpret results */

// 	ofstream ofs("output/test_tensor_matching.csv");
// 	ofs << "pair, score, hks_error1, hks_error2" << endl;
	vector<double> vTimes;
	for (int i = 0; i < 5; ++i) vTimes.push_back(10 * pow(1.4, i));
// 	for (int k = 0; k < vsize1; ++k)
// 	{
// 		int i1 = k;
// 		int i2 = (int)px[i1] - 1;
// 		int idx1 = vFeatures1[i1];
// 		int idx2 = vFeatures2[i2];
// 		ofs << idx1 << '-' << idx2 << ", " << pv[k] 
// 			 << ", " << calPointHKSSimilarity(pmp1, pmp2, idx1, idx2, vTimes, 0)
// 			 << ", " << calPointHKSSimilarity(pmp1, pmp2, idx1, idx2, vTimes, 1)
// 			 << endl;
// 	}	
// 	ofs.close();

	double result = ps[0];	// tensor matching score
//	result = 0.0;

	int count = 0;
	while(count++ < vsize1)
	{
		double *pmax = max_element(pv, pv+vsize1);	// max feature match score
		int imax = pmax - pv;	// shape 1 feature index
		MatchPair mpt;
		mpt.m_idx1 = vFeatures1[imax];
		int ind = (int)px[imax];		// shape 2 feature index matched to imax
		//if(ind<0 || ind>vsize2) continue;
		mpt.m_idx2 = vFeatures2[ind-1];
		mpt.m_score = *pmax;

		if (verbose)
		{
			cout << imax << "(" << mpt.m_idx1 << " - " << mpt.m_idx2 << "), score: " << *pmax
				 << ", Dissimilarity: " << calPointHksDissimilarity(pmp1, pmp2, mpt.m_idx1, mpt.m_idx2, vTimes, 1) 
				 << endl;
		}

		if(*pmax < para_thresh) break;

// 		if (calPointHksDissimilarity(pmp1, pmp2, mpt.m_idx1, mpt.m_idx2, vTimes, 1) >= 0.03) 
// 		{
// 			pv[imax] = 0.;
// 			continue;
// 		}

//		result += pv[imax];
		matched.push_back(mpt);

		pv[imax] = 0.0;
		// clear conflicted
		for(int i = 0; i < vsize1; i++)
		{
			if((int)px[i] == (int)px[imax])
				pv[i] = 0.0;
		}
	}

	mxDestroyArray(mxfeat1);
	mxDestroyArray(mxfeat2);
	mxDestroyArray(mxtris);
	mxDestroyArray(mxnumbs);
	// 	mxDestroyArray(mX2);
	// 	mxDestroyArray(vX2);
	// 	mxDestroyArray(score);

	return result;
}

double ShapeMatcher::TensorGraphMatching6( Engine *ep, 
        MeshHelper* pmp1, MeshHelper* pmp2, 
        const std::vector<HKSFeature>& vFeatures1, 
        const std::vector<HKSFeature>& vFeatures2, 
        std::vector<MatchPair>& matched, double para_t, 
        double para_thresh, bool verbose /*= false*/ )
{
	vector<int> vFeatIdx1(vFeatures1.size()), vFeatIdx2(vFeatures2.size());
	transform(vFeatures1.begin(), vFeatures1.end(), vFeatIdx1.begin(), [](const HKSFeature& feat){ return feat.m_index; });
	transform(vFeatures2.begin(), vFeatures2.end(), vFeatIdx2.begin(), [](const HKSFeature& feat){ return feat.m_index; });
	double matchScore = TensorGraphMatching6(ep, pmp1, pmp2, vFeatIdx1, vFeatIdx2, matched, para_t, para_thresh, verbose);
	return matchScore;
}	

double ShapeMatcher::TensorGraphMatching3(Engine *ep, 
    MeshHelper* pmp1, MeshHelper* pmp2, 
    const std::vector<int>& ct1, 
    const std::vector<int>& ct2, 
    std::vector<MatchPair>& matched, 
    double t, double thresh)
{
	// generate triangles
	int vsize1 = (int)ct1.size();	// input feature size 1
	int vsize2 = (int)ct2.size();  // input feature size 2

	vector<int> triangs;
	int i, j, k, tsize1, tsize2;

	// ***************************************
	// improve to local triangles
	// ***************************************/

	if(vsize1 > 8) 
	{
		for (i = 0; i < vsize1; i++)
		{
			for(j=0; j<8; j++)
			{
				if(i==j) continue;
				for(k=vsize1-1; k>vsize1-9; k--)
				{
					if(i==k || j==k) continue;
					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}
			}
		}
	}
	else
	{
		for(i=0; i<vsize1; i++)
		{
			for(j=0; j<vsize1; j++)
			{
				if(i==j) continue;
				for(k=0; k<vsize1; k++)
				{
					if(i==k || j==k) continue;
					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}
			}
		}
	}

	tsize1 = (int)triangs.size() / 3;
	tsize2 = vsize2 * vsize2 * vsize2;

	// compute feature descriptors

	mxArray *feat1, *feat2, *tris, *numbs, *mX2, *vX2, *score;
	feat1 = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *pfeat1 = mxGetPr(feat1);
	feat2 = mxCreateDoubleMatrix(3, tsize2, mxREAL);
	double *pfeat2 = mxGetPr(feat2);
	tris = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *ptris = mxGetPr(tris);
	numbs = mxCreateDoubleMatrix(1, 4, mxREAL);
	double *pnumbs = mxGetPr(numbs);

	pnumbs[0] = vsize1; // nP1
	pnumbs[1] = vsize2; // nP2
	pnumbs[2] = tsize1; // nT
	if(tsize1 > 60) pnumbs[3] = 60;     // nNN
	else pnumbs[3] = tsize1*0.5;

	for(i=0; i<tsize1*3; i++)
		ptris[i] = triangs[i];

	for(i=0; i<tsize1; i++)
	{
		int vi = triangs[i*3];
		int vj = triangs[i*3+1];
		int vk = triangs[i*3+2];
		ComputeTensorFeature3(pmp1, ct1[vi], ct1[vj], ct1[vk], t, &pfeat1[i*3]);
	}

	for(i=0; i<vsize2; i++)
	{
		for(j=0; j<vsize2; j++)
		{
			for(k=0; k<vsize2; k++)
			{
				ComputeTensorFeature3(pmp2, ct2[i], ct2[j], ct2[k], t, &pfeat2[((i*vsize2+j)*vsize2+k)*3]);
			}
		}
	}

	// invoke matlab for tensor matching

	engPutVariable(ep,"feat1",feat1);
	engPutVariable(ep,"feat2",feat2);
	engPutVariable(ep,"tris",tris);
	engPutVariable(ep,"numbs",numbs);
	//	system("PAUSE");	// for manually testing tensorMat
	engEvalString(ep, "[vX2,mX2,score]=tensorMat(feat1,feat2,tris,numbs);");
	vX2 = engGetVariable(ep, "vX2");
	double *pv = mxGetPr(vX2);	// best match of each shape-1 feature in the set of shape-2 features
	mX2 = engGetVariable(ep, "mX2");
	double *px = mxGetPr(mX2);  // corresponding match scores
	score = engGetVariable(ep, "score");
	double *ps = mxGetPr(score);

	// interpret results

	double result = ps[0];	// tensor matching score
	result = 0.0;

	int count = 0;
	while(count++ < vsize1)
	{
		double *pmax = max_element(pv, pv+vsize1);	// max feature match score
		int imax = pmax - pv;	// shape-1 feature index

		cout << imax << ": " << *pmax << endl;
		if(*pmax < thresh) 
			break;

		result += pv[imax];
		MatchPair mpt;
		mpt.m_idx1 = ct1[imax];
		int ind = (int)px[imax];		// shape-2 feature index matched to imax
		//if(ind<0 || ind>vsize2) continue;
		mpt.m_idx2 = ct2[ind-1];
		matched.push_back(mpt);

		pv[imax] = 0.0;
		// clear conflicted
		for(i=0; i<vsize1; i++)
		{
			if(px[i] == px[imax])
				pv[i] = 0.0;
		}
	}

	mxDestroyArray(feat1);
	mxDestroyArray(feat2);
	mxDestroyArray(tris);
	mxDestroyArray(numbs);
	// 	mxDestroyArray(mX2);
	// 	mxDestroyArray(vX2);
	// 	mxDestroyArray(score);

	return result;
}

void ShapeMatcher::readInRandPair( const std::string& filename )
{
	m_randPairs.clear();
	ifstream ifs(filename.c_str());

	double x, y;
	int count(0);
	ifs >> count;
	m_randPairs.reserve(count);
	for (int i = 0; i < count; ++i)
	{
		ifs >> x >> y;
		m_randPairs.push_back(make_pair(x, y));
		
	}
	cout << "## rand pair size: " << m_randPairs.size() << endl;
}

double ShapeMatcher::evaluateDistortion(const std::vector<MatchPair>& vIdMatchPair, CMesh* mesh1, CMesh* mesh2, const std::vector<std::pair<double, double> >& vRandPair, int rand_start /*= 0*/)
{
	int matchSize = vIdMatchPair.size(), rand_size = vRandPair.size();
	
	double distortSum = 0.;
	double avg_len_ratio = mesh1->getAvgEdgeLength() / mesh2->getAvgEdgeLength();

	const int total_run = 200;
	int count = 0;
	for (int k = rand_start; count < total_run || k < rand_size; ++k)
	{
		int pick1 = matchSize * vRandPair[k].first, pick2 = matchSize * vRandPair[k].second;
		int idx_11 = vIdMatchPair[pick1].m_idx1, idx_12 = vIdMatchPair[pick1].m_idx2;
		int idx_21 = vIdMatchPair[pick2].m_idx1, idx_22 = vIdMatchPair[pick2].m_idx2;

		if (idx_11 == idx_21) continue;

		double dist1 = calGeodesic(*mesh1, idx_11, idx_21), dist2 = calGeodesic(*mesh2, idx_12, idx_22) * avg_len_ratio;

		distortSum += abs(dist1 - dist2)/dist1;
		count++;
	}

	return distortSum / double(count);
}

void ShapeMatcher::forceInitialAnchors( const std::vector<MatchPair>& mp )
{
	m_nAlreadyMatchedLevel = m_nRegistrationLevels;
	m_vFeatureMatchingResults[m_nAlreadyMatchedLevel] = mp;
	m_bFeatureMatched = true;
}

void ShapeMatcher::evaluateWithGroundTruth( const std::vector<MatchPair>& vIdMatchPair )
{
	if (!m_bHasGroundTruth) {
		std::cout << "No ground truth available for evaluation!" << std::endl;
	}

	int count_valid = 0;
	int count_valid_strict = 0;
	for (auto iter = vIdMatchPair.begin(); iter != vIdMatchPair.end(); ++iter)
	{
		int idx1 = iter->m_idx1, idx2 = iter->m_idx2;
		if (mMatchGroundTruth.find(idx1) == mMatchGroundTruth.end()) continue;
		int idx1_map = mMatchGroundTruth[idx1];
		if (pOriginalMesh[1]->isInNeighborRing(idx1_map, idx2, 2))
			count_valid++;
		if (idx1_map == idx2)
			count_valid_strict++;
	}

	std::cout << "#Valid Matches: " << count_valid_strict << '/' << count_valid << '/' << vIdMatchPair.size() << endl;
}

void ShapeMatcher::evaluateWithGroundTruth(const MatchResult& result, const MatchResult& groundtruth, MatchEvaluation& eval) const
{
	const std::map<int,int>& resultMap = result.mMatchedPairs;
	const std::map<int,int>& groundTruthMap = groundtruth.mMatchedPairs;
	int countValid0 = 0, countValid1 = 0, countValid2 = 0;
	int largeErrorCount = 0;
	double errorSum(0);		
	int countCandidates(0);
	eval.mMatchedCount = resultMap.size();

	for (auto iter = resultMap.begin(); iter != groundTruthMap.end(); ++iter) {
		int idx1 = iter->first, idx2 = iter->second;
		if (groundTruthMap.find(idx1) == groundTruthMap.end()) continue;
		countCandidates++;

		int idxTrue = groundTruthMap.at(idx1);
		if (idx2 == idxTrue) { countValid0++; countValid1++; countValid2++; }
		else if (pOriginalMesh[1]->isInNeighborRing(idxTrue, idx2, 1)) { countValid1++; countValid2++; }
		else if (pOriginalMesh[1]->isInNeighborRing(idxTrue, idx2, 2)) { countValid2++; }

		double error = calGeodesic(*pOriginalMesh[1], idxTrue, idx2) / pOriginalMesh[1]->getAvgEdgeLength();
		if (error > 5) largeErrorCount++;
		errorSum += error;
	}
	errorSum /= countCandidates;

	eval.mMatch0Count = countValid0;
	eval.mMatch1Count = countValid1;
	eval.mMatch2Count = countValid2;
	eval.mLargeErrorCount = largeErrorCount;
	eval.mAvgErrInEdgeLength = errorSum;
}

const std::vector<MatchPair>& ShapeMatcher::getInitialMatchedFeaturePairs() const
{
	return m_vFeatureMatchingResults[m_nRegistrationLevels];
}

bool ShapeMatcher::loadInitialFeaturePairs( const std::string& filename )
{
	ifstream ifs(filename.c_str());
	if (!ifs) return false;
	int vsize;
	ifs >> vsize;
	vector<MatchPair> vPairs;
	for (int i = 0; i < vsize; ++i)
	{
		int idx1, idx2; double score;
		ifs >> idx1 >> idx2 >> score;
		vPairs.push_back(MatchPair(idx1, idx2, score));
	}

	forceInitialAnchors(vPairs);	
	return true;
}

void printVectorMatchPair(const vector<MatchPair>& vPairs, ostream& os)
{
	for_each(vPairs.begin(), vPairs.end(), [&](const MatchPair& mp) { os << '(' << mp.m_idx1 << ',' << mp.m_idx2 << ") ";});
	os << endl;
}

int countValidMatchPair(const vector<MatchPair>& vPairs)
{
	return count_if(vPairs.begin(), vPairs.end(), [](const MatchPair& mp){return mp.m_idx1 == mp.m_idx2;});
}

void ShapeMatcher::dataTesting1()
{
	vector<int> vFeatureID1, vFeatureID2;
	vFeatureID1.push_back(4297);
	vFeatureID1.push_back(4190);
	vFeatureID1.push_back(9527);
	vFeatureID1.push_back(6863);
	vFeatureID2.push_back(905);
	vFeatureID2.push_back(798);
	vFeatureID2.push_back(6135);
	vFeatureID2.push_back(4134);

	//vector<int> vSelected[2] = {tmesh1->getNeighborVertexIndex(5934, 1), tmesh2->getNeighborVertexIndex(5934, 1)};
	vector<int> vSelected[2] = {vFeatureID1, vFeatureID2};

	string filename = "output/register_test3.csv";
	ofstream ofs(filename.c_str());

	vector<double> vt;
	for (int i = 0; i < 20; ++i)
		vt.push_back(10 * pow(2.0, (double)i/2.));
	ofs << "Time";
	for (auto iter = vt.begin(); iter != vt.end(); ++iter)
		ofs << ", " << *iter;
	ofs << endl;

	for (int obj = 0; obj < 2; ++obj)
	{
		for (auto iter1 = vSelected[obj].begin(); iter1 != vSelected[obj].end(); ++iter1)
		{
			for (auto iter2 = vSelected[obj].begin(); iter2 != vSelected[obj].end(); ++iter2)
			{
				ofs << "Vertex" << obj+1 << "_" << *iter1 << '-' << *iter2;
				for (auto itert = vt.begin(); itert != vt.end(); ++itert)
				{
					double hk = calHK(pOriginalProcessor[obj], *iter1, *iter2, *itert);
					//hk = std::log(4.*PI*hk);
					//hk /= pOriginalProcessor[obj]->calHeatTrace(*iter2);
					ofs << ", " << hk;
				}
				ofs << endl;
			}
		}
	}	

	cout << "Collected data saved to \"" << filename << "\"" << endl;
}

double ShapeMatcher::calPointHksDissimilarity(MeshHelper* pmp1, MeshHelper* pmp2, int i1, int i2, const std::vector<double>& vTimes, int mode /*= 0*/)
{
	double errorSum(0);
	double maxError(0);
	for (auto iter = vTimes.begin(); iter != vTimes.end(); ++iter)
	{
		double hks1 = std::log(4*PI*(*iter) * calHK(pmp1, i1, i1, *iter));
		double hks2 = std::log(4*PI*(*iter) * calHK(pmp2, i2, i2, *iter));
		double error = abs(hks1-hks2);
		errorSum += error*error;
		maxError = max(maxError, error);
	}
	if (mode == 0) return maxError;
	else return sqrt(errorSum / vTimes.size());
}



void ShapeMatcher::SimplePointMatching(MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, const std::vector<double>& vTimes, std::vector<MatchPair>& matchedResult, bool verbose /*= false*/ )
{
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);
	matchedResult.clear();

	int vsize1 = vFeatures1.size(), vsize2 = vFeatures2.size();
	int tsize = vTimes.size();
	vector<VecNd> vSig1(vsize1), vSig2(vsize2);

	vector<double> vTrace1(tsize), vTrace2(tsize);
	for (int s = 0; s < tsize; ++s)
	{
        vTrace1[s] = calHeatTrace(es1, vTimes[s]);
        vTrace2[s] = calHeatTrace(es2, vTimes[s]);
	}

	for (int i = 0; i < vsize1; ++i) 
	{
		vSig1[i].resize(tsize);
		for (int j = 0; j < tsize; ++j)
			vSig1[i][j] = calHK(pmp1, vFeatures1[i], vFeatures1[i], vTimes[j]) / vTrace1[j];
	}
	for (int i = 0; i < vsize2; ++i) 
	{
		vSig2[i].resize(tsize);
		for (int j = 0; j < tsize; ++j)
			vSig2[i][j] = calHK(pmp2, vFeatures2[i], vFeatures2[i], vTimes[j]) / vTrace2[j];
	}
	
	double *hksSim = new double[vsize1 * vsize2];

	for (int i = 0; i < vsize1; ++i)
	{
		for (int j = 0; j <= i; ++j)
			hksSim[i * vsize2 + j] = std::pow((vSig1[i] - vSig2[j]).norm2(), 2);
		for (int j = i+1; j < vsize2; ++j)
			hksSim[i* vsize2 + j] = 1e10;
	}

	double vMin = 1e10 - 1;
	while (vMin < 1e10 - 0.1)
	{
		double *pMin = min_element(hksSim, hksSim + vsize1 * vsize2);
		vMin = *pMin;
		size_t pos = distance(hksSim, pMin);
		int i2 = pos % vsize2;
		int i1 = pos / vsize2;
		matchedResult.push_back(MatchPair(vFeatures1[i1], vFeatures2[i2], vMin));

		if (verbose)
		{
			cout << matchedResult.size() << "(" << matchedResult.back().m_idx1 << " - " << matchedResult.back().m_idx2 << "), score: " << matchedResult.back().m_score	<< endl;
		}

		for (int k = 0; k < vsize1; ++k)
			hksSim[k*vsize2 + i2] = 1e10;
		for (int k = 0; k < vsize2; ++k)
			hksSim[i1*vsize2 + k] = 1e10;
	}
	
	delete []hksSim;
}

void ShapeMatcher::matchFeatureSimple()
{
	const CMesh *mesh1 = pOriginalMesh[0], *mesh2 = pOriginalMesh[1];
	const vector<HKSFeature>& vftFine1 = m_vFeatures[0];
	const vector<HKSFeature>& vftFine2 = m_vFeatures[1];
	vector<MatchPair> vPairs;
	vector<int> vFeatures1(vftFine1.size()), vFeatures2(vftFine2.size());
	transform(vftFine1.begin(), vftFine1.end(), vFeatures1.begin(), [](const HKSFeature& feat){ return feat.m_index; });
	transform(vftFine2.begin(), vftFine2.end(), vFeatures2.begin(), [](const HKSFeature& feat){ return feat.m_index; });
	
	vector<double> vTimes;
	for(int i = 0; i < 10; ++i)
		vTimes.push_back(10. * pow(2., i/2.));

	SimplePointMatching(pOriginalProcessor[0], pOriginalProcessor[1], vFeatures1, vFeatures2, vTimes, vPairs, true);

	forceInitialAnchors(vPairs);
}

double ShapeMatcher::TensorMatchingExt( Engine *ep, MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPairs, int highOrderFeatureType, double vPara[], std::ostream& logout, bool verbose /*= false*/ )
{
	/******************************
	* featureType = 0 ---- original angle based triangle construction
	** vPara[0]: timescale
	** vPara[1]: thresh

	* featureType = 1 ---- 6-tuple, 3 hks + 3 hk
	** vPara[0]: timescale
	** vPara[1]: thresh

	* featureType = 2 ---- 6-tuple, 3 hk + 3 hk to anchor
	** vPara[0]: timescale
	** vPara[1]: thresh
	** vPara[2]: anchor 1-1
	** vPara[3]: anchor 1-2

	* featureType = 3 ---- 12-tuple
	** vPara[0]: timescale 1
	** vPara[1]: thresh
	** vPara[2]: timescale 2

	******************************/
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);

	vMatchedPairs.clear();
	const int vsize1 = (int)vFeatures1.size();	// input feature size 1
	const int vsize2 = (int)vFeatures2.size();    // input feature size 2
	if (vsize1 < 3 || vsize2 < 3) return 0;

	if (vsize1 > vsize2)
	{
		vector<MatchPair> oppositeMatched;
		double score;
		if (highOrderFeatureType == 2)
		{
			double vParaOpp[] = {vPara[0], vPara[1], vPara[3], vPara[2]};
			score = TensorMatchingExt(ep, pmp2, pmp1, vFeatures2, vFeatures1, oppositeMatched, highOrderFeatureType, vParaOpp, logout);
		}
		else 
			score = TensorMatchingExt(ep, pmp2, pmp1, vFeatures2, vFeatures1, oppositeMatched, highOrderFeatureType, vPara, logout);

		vMatchedPairs.resize(oppositeMatched.size());
		transform(oppositeMatched.begin(), oppositeMatched.end(), vMatchedPairs.begin(), [](const MatchPair& mp){return MatchPair(mp.m_idx2, mp.m_idx1, mp.m_score);});
		return score;
	}

	// generate triangles
	vector<int> triangs;
	if (highOrderFeatureType == -1)
	{
// 		if(vsize1>8) 
// 		{
// 			for(int i=0; i<vsize1; i++)
// 			{
// 				for(int j=0; j<8; j++)
// 				{
// 					if(i==j) continue;
// 					for(int k=vsize1-1; k>vsize1-9; k--)
// 					{
// 						if(i==k || j==k) continue;
// 						triangs.push_back(i);
// 						triangs.push_back(j);
// 						triangs.push_back(k);
// 					}
// 				}
// 			}
// 		}
// 		else
// 		{
// 			for(int i=0; i<vsize1; i++)
// 			{
// 				for(int j=0; j<vsize1; j++)
// 				{
// 					if(i==j) continue;
// 					for(int k=0; k<vsize1; k++)
// 					{
// 						if(i==k || j==k) continue;
// 						triangs.push_back(i);
// 						triangs.push_back(j);
// 						triangs.push_back(k);
// 					}
// 				}
// 			}
// 		}
		for (int i = 0; i < vsize1; i++)
		{
			for (int j = 0; j < vsize1; j++)
			{
				if (i == j) continue;
				for (int k = 0; k < vsize1; k++)
				{
					if (i == k || k == j) continue;
					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}			
			}
		}
	}
	else if (highOrderFeatureType == 1 || highOrderFeatureType == 2 || highOrderFeatureType == 0)
	{
		double t = vPara[0];
		for (int i = 0; i < vsize1; i++)
		{
			for (int j = i+1; j < vsize1; j++)
			{
				double c1 = calHK(pmp1, vFeatures1[i], vFeatures1[j], t);
				for (int k = j+1; k < vsize1; k++)
				{
					double c2 = calHK(pmp1, vFeatures1[j], vFeatures1[k], t),
						   c3 = calHK(pmp1, vFeatures1[i], vFeatures1[k], t);
					int count_valid = 0;
					if (c1 >= 1e-5) count_valid++;
					if (c2 >= 1e-5) count_valid++;
					if (c3 >= 1e-5) count_valid++;

					if (count_valid == 0) continue;

					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}			
			}
		}
	}
	else if (highOrderFeatureType == 3)
	{
		double t1 = vPara[0], t2 = vPara[2];
		for (int i = 0; i < vsize1; i++)
		{
			for (int j = i+1; j < vsize1; j++)
			{
				double c1 = calHK(pmp1, vFeatures1[i], vFeatures1[j], t1);
				double d1 = calHK(pmp1, vFeatures1[i], vFeatures1[j], t2);
				for (int k = j+1; k < vsize1; k++)
				{
					double c2 = calHK(pmp1, vFeatures1[j], vFeatures1[k], t1),
						   c3 = calHK(pmp1, vFeatures1[i], vFeatures1[k], t1);
					double d2 = calHK(pmp1, vFeatures1[j], vFeatures1[k], t2),
						   d3 = calHK(pmp1, vFeatures1[i], vFeatures1[k], t2);
					int count_valid1 = 0;
					if (c1 >= 1e-5) count_valid1++;
					if (c2 >= 1e-5) count_valid1++;
					if (c3 >= 1e-5) count_valid1++;
					int count_valid2 = 0;
					if (d1 >= 1e-5) count_valid2++;
					if (d2 >= 1e-5) count_valid2++;
					if (d3 >= 1e-5) count_valid2++;

					if (count_valid1 == 0 || count_valid2 == 0) continue;

					triangs.push_back(i);
					triangs.push_back(j);
					triangs.push_back(k);
				}			
			}
		}

		// 	for (int i = 0; i < vsize1; i++)
		// 	{
		// 		for (int j = 0; j < vsize1; j++)
		// 		{
		// 			for (int k = 0; k < vsize1; k++)
		// 			{
		// 				if (i == j || i == k || k == j) continue;
		// 				triangs.push_back(i);
		// 				triangs.push_back(j);
		// 				triangs.push_back(k);
		// 			}			
		// 		}
		// 	}
	}

	if (triangs.size() < 3)
	{
		logout << "no valid triangles!" << endl;
		return 0;
	}

	int tsize1 = (int)triangs.size() / 3;

	if (verbose)
		logout << "#Total query triangles: " << tsize1 << endl;
	
	int tsize2 = vsize2 * vsize2 * vsize2;

	// compute feature descriptors
	int FeatureDim(3);
	switch(highOrderFeatureType)
	{
	case 0: FeatureDim = 3; break;
	case 1: FeatureDim = 6; break;
	case 2: FeatureDim = 6; break;
	case 3: FeatureDim = 12; break;
	}
	
	mxArray *mxfeat1, *mxfeat2, *mxtris, *mxnumbs, *mX2, *vX2, *score;
	mxfeat1 = mxCreateDoubleMatrix(FeatureDim, tsize1, mxREAL);
	double *pfeat1 = mxGetPr(mxfeat1);
	mxfeat2 = mxCreateDoubleMatrix(FeatureDim, tsize2, mxREAL);
	double *pfeat2 = mxGetPr(mxfeat2);
	mxtris = mxCreateDoubleMatrix(3, tsize1, mxREAL);
	double *ptris = mxGetPr(mxtris);
	mxnumbs = mxCreateDoubleMatrix(1, 4, mxREAL);
	double *pnumbs = mxGetPr(mxnumbs);

	const int maxNN = 32;
	pnumbs[0] = vsize1; // nP1
	pnumbs[1] = vsize2; // nP2
	pnumbs[2] = tsize1; // nT
	if(tsize1 > 2*maxNN) pnumbs[3] = maxNN;     // nNN
	else pnumbs[3] = max(1, tsize1 / 2);

	copy(triangs.begin(), triangs.end(), ptris);


	//	for(int i = 0; i < tsize1; i++)
	switch(highOrderFeatureType)
	{
	case 0:
		{
			Concurrency::parallel_for(0, tsize1, [&](int i)
			{
				int vi = triangs[i*3];
				int vj = triangs[i*3+1];
				int vk = triangs[i*3+2];
				ComputeTensorFeature3(pmp1, vFeatures1[vi], vFeatures1[vj], vFeatures1[vk], vPara[0], pfeat1 + i*FeatureDim);
			}
			);
			Concurrency::parallel_for(0, vsize2, [&](int i)
			{
				for(int j = 0; j < vsize2; j++)
				{
					for(int k = 0; k < vsize2; k++)
					{
						ComputeTensorFeature3(pmp2, vFeatures2[i], vFeatures2[j], vFeatures2[k], vPara[0], pfeat2 + ((i*vsize2+j)*vsize2+k)*FeatureDim);
					}
				}
			}
			);
			break;
		}
	case 1:
		{
			Concurrency::parallel_for(0, tsize1, [&](int i)
			{
				int vi = triangs[i*3];
				int vj = triangs[i*3+1];
				int vk = triangs[i*3+2];
				ComputeTensorFeature6(pmp1, vFeatures1[vi], vFeatures1[vj], vFeatures1[vk], vPara[0], pfeat1 + i*FeatureDim);
			}
			);
			Concurrency::parallel_for(0, vsize2, [&](int i)
			{
				for(int j = 0; j < vsize2; j++)
				{
					for(int k = 0; k < vsize2; k++)
					{
						ComputeTensorFeature6(pmp2, vFeatures2[i], vFeatures2[j], vFeatures2[k], vPara[0], pfeat2 + ((i*vsize2+j)*vsize2+k)*FeatureDim);
					}
				}
			}
			);
			break;
		}
	case 2:
		{
			Concurrency::parallel_for(0, tsize1, [&](int i)
			{
				int vi = triangs[i*3];
				int vj = triangs[i*3+1];
				int vk = triangs[i*3+2];
				ComputeTensorFeatureAnchor(pmp1, vFeatures1[vi], vFeatures1[vj], vFeatures1[vk], vPara[2], vPara[0], pfeat1 + i*FeatureDim);
			}
			);
			Concurrency::parallel_for(0, vsize2, [&](int i)
			{
				for(int j = 0; j < vsize2; j++)
				{
					for(int k = 0; k < vsize2; k++)
					{
						ComputeTensorFeatureAnchor(pmp2, vFeatures2[i], vFeatures2[j], vFeatures2[k], vPara[3], vPara[0], pfeat2 + ((i*vsize2+j)*vsize2+k)*FeatureDim);
					}
				}
			}
			);
			break;
		}
	case 3:
		{
			Concurrency::parallel_for(0, tsize1, [&](int i)
			{
				int vi = triangs[i*3];
				int vj = triangs[i*3+1];
				int vk = triangs[i*3+2];
				ComputeTensorFeature12(pmp1, vFeatures1[vi], vFeatures1[vj], vFeatures1[vk], vPara[0], vPara[2], pfeat1 + i*FeatureDim);
			}
			);
			Concurrency::parallel_for(0, vsize2, [&](int i)
			{
				for(int j = 0; j < vsize2; j++)
				{
					for(int k = 0; k < vsize2; k++)
					{
						ComputeTensorFeature12(pmp2, vFeatures2[i], vFeatures2[j], vFeatures2[k], vPara[0], vPara[2], pfeat2 + ((i*vsize2+j)*vsize2+k)*FeatureDim);
					}
				}
			}
			);
			break;
		}
	}

	// invoke matlab for tensor matching

#define HK_NORMALIZE
#ifdef HK_NORMALIZE
	if (highOrderFeatureType == 1 || highOrderFeatureType == 2)	// scaling hk on different models
	{
		double ht1 = calHeatTrace(es1, vPara[0]), ht2 = calHeatTrace(es2, vPara[0]);
		transform(pfeat1, pfeat1+FeatureDim*tsize1, pfeat1, [=](double v){ return v/ht1; });
		transform(pfeat2, pfeat2+FeatureDim*tsize2, pfeat2, [=](double v){ return v/ht2; });
	}
#endif	

	engPutVariable(ep, "feat1", mxfeat1);
	engPutVariable(ep, "feat2", mxfeat2);
	engPutVariable(ep, "tris", mxtris);
	engPutVariable(ep, "numbs", mxnumbs);
	//	system("PAUSE");	// for manually testing tensorMat
	engEvalString(ep, "[vX2,mX2,score]=tensorMat(feat1,feat2,tris,numbs);");
	mX2 = engGetVariable(ep, "mX2");
	double *px = mxGetPr(mX2);  // best match of each shape-1 feature in the set of shape-2 features
	vX2 = engGetVariable(ep, "vX2");
	double *pv = mxGetPr(vX2);	// corresponding match scores
	score = engGetVariable(ep, "score");	// general matching score
	double *ps = mxGetPr(score);

	/* interpret results */

	// 	ofstream ofs("output/test_tensor_matching.csv");
	// 	ofs << "pair, score, hks_error1, hks_error2" << endl;
	vector<double> vTimes;
	for (int i = 0; i < 5; ++i) vTimes.push_back(10 * pow(1.4, i));
	// 	for (int k = 0; k < vsize1; ++k)
	// 	{
	// 		int i1 = k;
	// 		int i2 = (int)px[i1] - 1;
	// 		int idx1 = vFeatures1[i1];
	// 		int idx2 = vFeatures2[i2];
	// 		ofs << idx1 << '-' << idx2 << ", " << pv[k] 
	// 			 << ", " << calPointHKSSimilarity(pmp1, pmp2, idx1, idx2, vTimes, 0)
	// 			 << ", " << calPointHKSSimilarity(pmp1, pmp2, idx1, idx2, vTimes, 1)
	// 			 << endl;
	// 	}	
	// 	ofs.close();

	double result = ps[0];	// tensor matching score
	//	result = 0.0;

	int count = 0;
	while(count++ < vsize1)
	{
		double *pmax = max_element(pv, pv+vsize1);	// max feature match score
		int imax = pmax - pv;	// shape 1 feature index
		MatchPair mpt;
		mpt.m_idx1 = vFeatures1[imax];
		int ind = (int)px[imax];		// shape 2 feature index matched to imax
		//if(ind<0 || ind>vsize2) continue;
		mpt.m_idx2 = vFeatures2[ind-1];
		mpt.m_score = *pmax;

		if (verbose)
			logout << imax << "(" << mpt.m_idx1 << " - " << mpt.m_idx2 << "), score: " << *pmax
			//	<< ", Dissimilarity: " << calPointHksDissimilarity(pmp1, pmp2, mpt.m_idx1, mpt.m_idx2, vTimes, 1) 
				<< endl;
		

		if(*pmax < vPara[1]) break;

		// 		if (calPointHksDissimilarity(pmp1, pmp2, mpt.m_idx1, mpt.m_idx2, vTimes, 1) >= 0.03) 
		// 		{
		// 			pv[imax] = 0.;
		// 			continue;
		// 		}

		//		result += pv[imax];
		vMatchedPairs.push_back(mpt);

		pv[imax] = 0.0;
		// clear conflicted
		for(int i = 0; i < vsize1; i++)
		{
			if((int)px[i] == (int)px[imax])
				pv[i] = 0.0;
		}
	}

	mxDestroyArray(mxfeat1);
	mxDestroyArray(mxfeat2);
	mxDestroyArray(mxtris);
	mxDestroyArray(mxnumbs);
	// 	mxDestroyArray(mX2);
	// 	mxDestroyArray(vX2);
	// 	mxDestroyArray(score);

	return result;
}

void ShapeMatcher::ComputeTensorFeatureAnchor(MeshHelper* pmp, int i, int j, int k, int origin, double t, double* sang )
{
	double d1 = calHK(pmp, i, j, t);
	double d2 = calHK(pmp, j, k, t);
	double d3 = calHK(pmp, k, i, t);
	double s1 = calHK(pmp, origin, i, t);
	double s2 = calHK(pmp, origin, j, t);
	double s3 = calHK(pmp, origin, k, t);

	sang[0] = d1;
	sang[1] = d2;
	sang[2] = d3;
	sang[3] = s1;
	sang[4] = s2;
	sang[5] = s3;
}

void ShapeMatcher::ComputeTensorFeature12(MeshHelper* pmp, int i, int j, int k, double t1, double t2, double* sang )
{
    sang[0] = calHK(pmp, i, j, t1);
    sang[1] = calHK(pmp, j, k, t1);
    sang[2] = calHK(pmp, k, i, t1);
    sang[3] = calHK(pmp, i, i, t1);
    sang[4] = calHK(pmp, j, j, t1);
    sang[5] = calHK(pmp, k, k, t1);
    sang[6] = calHK(pmp, i, j, t2);
	sang[7] = calHK(pmp, j, k, t2);
	sang[8] = calHK(pmp, k, i, t2);
	sang[9] = calHK(pmp, i, i, t2);
	sang[10] = calHK(pmp, j, j, t2);
	sang[11] = calHK(pmp, k, k, t2);
}

double ShapeMatcher::PairMatchingExt(Engine* ep, MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, int method, double vPara[], std::ostream& logout, bool verbose /*= false*/ )
{
	/**** method = 0
	 ** vPara[0] = thresh1
	 ** vPara[1] = thresh2
	 ** vPara[2] = thresh3
	*/
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);
	vMatchedPair.clear();
	int size1 = (int)vFeatures1.size();
	int size2 = (int)vFeatures2.size();
	double thresh1 = vPara[0];	//0.52
	double thresh2 = vPara[1];  //0.02, 0.1
	double thresh3 = vPara[2];  //0.01, 0.1

	vector<MatchPair> vTmpMatchPairs;
	
	double sigma1 = 4.0 * thresh1;
	double sigma2 = thresh2; //0.02;	//0.02; 0.1

	vector<double> vTimes;
	vector<double> vHT1, vHT2;
	const double tl = 10, tn = 8;
	for (int i = 0; i < tn; ++i)
	{
		vTimes.push_back(tl * pow(2.0,i));
        vHT1.push_back(calHeatTrace(es1, vTimes.back()));
		vHT2.push_back(calHeatTrace(es2, vTimes.back()));
	}

	double *tmpScores = new double[size1*size2];
//	for (int i = 0; i < size1; ++i)
	Concurrency::parallel_for(0, size1, [&](int i)
	{
		for(int j = 0; j < size2; j++)
		{
			VecNd v1(tn), v2(tn);
			for(int k = 0; k < tn; k++)
			{
				v1[k] = calHK(pmp1, vFeatures1[i], vFeatures1[i], vTimes[k]) / vHT1[k];
				v2[k] = calHK(pmp2, vFeatures2[j], vFeatures2[j], vTimes[k]) / vHT2[k];
			}
			double dist = (v1 - v2).norm2();
			double d = dist * dist / tn;
			tmpScores[i*size2+j] = d;
		}
	}
	);

	for (int i = 0; i < size1; ++i)
	{
		for (int j = 0; j < size2; ++j)
		{
			if (tmpScores[i*size2+j] < thresh1)
			{
				double score = std::exp(-tmpScores[i*size2+j]/sigma1);
				vTmpMatchPairs.push_back(MatchPair(vFeatures1[i], vFeatures2[j], tl, tn, score));
			}
		}
	}
	delete []tmpScores;

	/* ---- create affinity matrix (compatibility of each candidate match) ---- */
	const int affinitySize = (int)vTmpMatchPairs.size();
	mxArray *AM, *VM, *VA;
	AM = mxCreateDoubleMatrix(affinitySize, affinitySize, mxREAL);
	double *am = mxGetPr(AM);
	Concurrency::parallel_for(0, affinitySize, [&](int i)
	{
		am[i*affinitySize+i] = vTmpMatchPairs[i].m_score; // diagonal

		for(int j = i+1; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i].m_idx1 == vTmpMatchPairs[j].m_idx1)^(vTmpMatchPairs[i].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor conflict; 1-to-multiple not desired? 
			{
				am[i*affinitySize+j] = am[j*affinitySize+i] = 0.0;
				continue;
			}

			int x1 = vTmpMatchPairs[i].m_idx1;
			int x2 = vTmpMatchPairs[i].m_idx2;
			int y1 = vTmpMatchPairs[j].m_idx1;
			int y2 = vTmpMatchPairs[j].m_idx2;

			VecNd v1(tn), v2(tn);
			double t = tl;
			for(int k = 0; k < tn; k++)
			{
				v1[k] = calHK(pmp1, x1, y1, vTimes[k]) / vHT1[k];
				v2[k] = calHK(pmp2, x2, y2, vTimes[k]) / vHT2[k];
			}
			double dist = (v1 - v2).norm2();
			double ds = dist * dist / tn;
			am[i*affinitySize+j] = am[j*affinitySize+i] = exp(-ds/sigma2);
		}
	});

	/* ---- solving the greatest eigenvector (PCA?) ---- */
	engPutVariable(ep, "AM", AM);
	engEvalString(ep, "[VM,VA] = spectral_embedding(AM);");	//computing leading eigenvector of A using svd
	VM = engGetVariable(ep, "VM");
	double *vm = mxGetPr(VM);
	VA = engGetVariable(ep, "VA");
	double *va = mxGetPr(VA);		//?? not referenced?

	const double c_thresh = thresh3;	//0.01;	//was 0.1
	std::vector<MatchPair> mpc1;

	std::transform(vm, vm+affinitySize, vm, [](double v){ return abs(v); });

	double totalMatchScore(0.);

	while(1)
	{
		int i_max = -1;
		double v_max = -1;
		for(int i = 0; i < affinitySize; i++)
		{
			if(vm[i] > v_max)
			{
				v_max = vm[i];
				i_max = i;
			}
		}

		if (v_max <= c_thresh) 
		{
			if (verbose) logout << "-- Discarded v_max: " << v_max << endl;
			break;
		}

		totalMatchScore += v_max;

		const int curMatchSize = (int)mpc1.size();		
		if (verbose) logout << "\t" << curMatchSize << ": " << i_max << ',' << v_max << endl;

		vm[i_max] = 0.0;
		mpc1.push_back(vTmpMatchPairs[i_max]);

		// now a max candidate is found, remove conflicting candidates
		for(int j = 0; j < affinitySize; j++)
		{
			if((vTmpMatchPairs[i_max].m_idx1 == vTmpMatchPairs[j].m_idx1) ^ (vTmpMatchPairs[i_max].m_idx2 == vTmpMatchPairs[j].m_idx2)) // xor, conflict
				vm[j] = 0.0;
		}
	} // end of while()

	mxDestroyArray(VM);
	mxDestroyArray(VA);
	mxDestroyArray(AM);

	vMatchedPair = mpc1;
	if (verbose) logout << "Pair graph matching computed!" << endl;

	return totalMatchScore;
}

void ShapeMatcher::HKSMatchingExt(MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, int method, double vPara[], std::ostream& logout, bool verbose /*= false*/ )
{
	/*********************************
	1. method = 0 (HKS similarity)
	** vPara[0] = tl
	** vPara[1] = tu
	** vPara[2] = thresh
	2. method = 1 (HKM)
	** vPara[0] = tl
	** vPara[1] = tu
	** vPara[2] = thresh
	** vPara[3] = anchor1
	** vPara[4] = anchor2 (may not be necessary)
	*********************************/
    EigenSystem &es1 = pmp1->getEigenSystem(CotFormula), &es2 = pmp2->getEigenSystem(CotFormula);
	double tl = vPara[0];
	double tu = vPara[1];
	double thresh = vPara[2];	//0.1

	vMatchedPair.clear();

	if (tu < tl) return;
	int tn = int(std::log(tu/tl)/std::log(2.0)) + 1;	//only consider the overlapping times
	vector<double> vTimes;
	for (int i = 0; i < tn; ++i)
	{
		vTimes.push_back(tl * pow(2.,i));
	}
	
	int vsize1 = vFeatures1.size(), vsize2 = vFeatures2.size();
	int tsize = vTimes.size();
	vector<VecNd> vSig1(vsize1), vSig2(vsize2);

	vector<double> vTrace1(tsize), vTrace2(tsize);
	for (int s = 0; s < tsize; ++s)
	{
        vTrace1[s] = calHeatTrace(es1, vTimes[s]);
        vTrace2[s] = calHeatTrace(es2, vTimes[s]);
	}

	if (method == 0)
	{
		for (int i = 0; i < vsize1; ++i) 
		{
			vSig1[i].resize(tsize);
			for (int j = 0; j < tsize; ++j)
				vSig1[i][j] = calHK(pmp1, vFeatures1[i], vFeatures1[i], vTimes[j]) / vTrace1[j];
		}
		for (int i = 0; i < vsize2; ++i) 
		{
			vSig2[i].resize(tsize);
			for (int j = 0; j < tsize; ++j)
				vSig2[i][j] = calHK(pmp2, vFeatures2[i], vFeatures2[i], vTimes[j]) / vTrace2[j];
		}
	}
	else if (method == 1)
	{
		int anchor1 = (int)vPara[3];
		int anchor2 = (int)vPara[4];
		for (int i = 0; i < vsize1; ++i) 
		{
			vSig1[i].resize(tsize);
			for (int j = 0; j < tsize; ++j)
				vSig1[i][j] = calHK(pmp1, anchor1, vFeatures1[i], vTimes[j]) /*/ vTrace1[j]*/;
		}
		for (int i = 0; i < vsize2; ++i) 
		{
			vSig2[i].resize(tsize);
			for (int j = 0; j < tsize; ++j)
				vSig2[i][j] = calHK(pmp2, anchor2, vFeatures2[i], vTimes[j]) /*/ vTrace2[j]*/;
		}
	}

	double *hksSim = new double[vsize1 * vsize2];

	for (int i = 0; i < vsize1; ++i)
	{
		for (int j = 0; j < vsize2; ++j)
			hksSim[i * vsize2 + j] = vSig1[i].distEculidean2(vSig2[j]) / tsize;
	}

	double vMin = 1e10 - 1;
	while (vMin < 1e10 - 0.1)
	{
		double *pMin = min_element(hksSim, hksSim + vsize1 * vsize2);
		vMin = *pMin;

		if (vMin > thresh) break;

		size_t pos = pMin - hksSim;
		int i2 = pos % vsize2;
		int i1 = pos / vsize2;
		vMatchedPair.push_back(MatchPair(vFeatures1[i1], vFeatures2[i2], vMin));

		if (verbose)
		{
			logout << vMatchedPair.size() << "(" << vMatchedPair.back().m_idx1 << " - " << vMatchedPair.back().m_idx2 << "), score: " << vMatchedPair.back().m_score	<< endl;
		}

		for (int k = 0; k < vsize1; ++k)
			hksSim[k*vsize2 + i2] = 1e10;
		for (int k = 0; k < vsize2; ++k)
			hksSim[i1*vsize2 + k] = 1e10;
	}

	delete []hksSim;
}

void ShapeMatcher::sparseMatchingTesting()
{
	int vsize = pOriginalMesh[0]->vertCount();
	double outlier_ratio = 0.2;
	ofstream ofs("output/evaluate_hkt.txt");
	ofstream ofs2("output/result_hkt.csv");

	vector<int> vN;
	for (int n = 40; n <= 75; n+=5) vN.push_back(n);
	const int groupNum = 8;
	for (auto iter_N = vN.begin(); iter_N != vN.end(); ++iter_N)
	{
		cout << "## N = " << *iter_N << " ##" << endl;
		ofs << "## N = " << *iter_N << " ##" << endl;
		const int N = *iter_N;

		/*   generate random test input    */
		vector<vector<int> > vvInput1(groupNum), vvInput2(groupNum);
		for (int group = 0; group < groupNum; ++group)
		{
			ofs << "  group #" << group << endl;

			vector<int> vToMatch(N), vOutlier1(outlier_ratio*N), vOutlier2(outlier_ratio*N);
			generate(vToMatch.begin(), vToMatch.end(), [&](){ return double(rand())/double(RAND_MAX)*vsize;});
			generate(vOutlier1.begin(), vOutlier1.end(), [&](){ return double(rand())/double(RAND_MAX)*vsize;});
			generate(vOutlier2.begin(), vOutlier2.end(), [&](){ return double(rand())/double(RAND_MAX)*vsize;});

			vvInput1[group].insert(vvInput1[group].begin(), vToMatch.begin(), vToMatch.end());
			vvInput2[group].insert(vvInput2[group].begin(), vToMatch.begin(), vToMatch.end());
			vvInput1[group].insert(vvInput1[group].end(), vOutlier1.begin(), vOutlier1.end());
			vvInput2[group].insert(vvInput2[group].end(), vOutlier2.begin(), vOutlier2.end());

			ofs << "  input " << group << "-1: ";
			for_each(vvInput1[group].begin(), vvInput1[group].end(), [&](int v){ ofs << v << ' '; });
			ofs << '\n';
			ofs << "  input " << group << "-2: ";
			for_each(vvInput2[group].begin(), vvInput2[group].end(), [&](int v){ ofs << v << ' '; });
			ofs << '\n';
		}
		
		double sum_ratio_matched[8] = {0.}, sum_ratio_correct[8]= {0.};
		int test_cat = 0;
		vector<MatchPair> vMatched;
		int validMatched;
		CStopWatch timer;

		/*-------------------------------------------------------------*/
		bool gotya = false;
		timer.startTimer();
		int selectedGroup = 0;
		for (int group = 0; group < groupNum; ++group)
		{
			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
			double vPara[] = {40, 0.7, 80};
			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 3, vPara, cout);
			ofs << "  Match results(method=3, t=40,80): "; 
			printVectorMatchPair(vMatched, ofs);
			validMatched = countValidMatchPair(vMatched);

			if ((double)validMatched / (double)N > 0.9)
			{
				for (auto iter = vMatched.begin(); iter != vMatched.end(); ++iter)
				{
					if (iter->m_idx1 == iter->m_idx2) iter->m_note = 1;
					else iter->m_note = -1;
				}
				forceInitialAnchors(vMatched);

				MeshFeatureList* mfl1 =new MeshFeatureList;
				MeshFeatureList* mfl2 = new MeshFeatureList;
				for (int i = 0; i < N * 6 / 5; ++i)
				{
					MeshFeature* f1 = new MeshFeature;
					f1->m_index = vinput1[i];
					if (i < N) f1->m_note = 1;
					else f1->m_note = -1;
					mfl1->addFeature(f1);

					MeshFeature* f2 = new MeshFeature;
					f2->m_index = vinput2[i];
					if (i < N) f2->m_note = 1;
					else f2->m_note = -1;
					mfl2->addFeature(f2);
				}				
				
				pOriginalProcessor[0]->getMesh()->addAttrMeshFeatures(*mfl1, StrAttrFeatureUnnamed);
				pOriginalProcessor[1]->getMesh()->addAttrMeshFeatures(*mfl2, StrAttrFeatureUnnamed);


				selectedGroup = group;

				gotya = true;
				break;
			}

			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
		}
		timer.stopTimer();
		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
		test_cat++;

		if (gotya)
		{
			{
				int group = selectedGroup;
				const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
				int anchor = 1210;//double(rand())/double(RAND_MAX) * vsize;
				double vPara[] = {10, 1280, 0.52, anchor};
				HKSMatchingExt(pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 1, vPara, cout, false);
				
				MeshFeatureList* mfl1 =new MeshFeatureList;
				MeshFeatureList* mfl2 = new MeshFeatureList;
				for (int i = 0; i < N * 6 / 5; ++i)
				{
					MeshFeature* f1 = new MeshFeature;
					f1->m_index = vinput1[i];
					if (i < N) f1->m_note = 1;
					else f1->m_note = -1;
					mfl1->addFeature(f1);

					MeshFeature* f2 = new MeshFeature;
					f2->m_index = vinput2[i];
					if (i < N) f2->m_note = 1;
					else f2->m_note = -1;
					mfl2->addFeature(f2);
				}				

				MeshFeature* anchor1 = new MeshFeature;
				anchor1->m_index = anchor;
				anchor1->m_note = 0;
				mfl1->addFeature(anchor1);
				MeshFeature* anchor2 = new MeshFeature;
				anchor2->m_index = anchor;
				anchor2->m_note = 0;
				mfl2->addFeature(anchor2);
				pOriginalProcessor[0]->getMesh()->addAttrMeshFeatures(*mfl1, StrAttrFeatureUnnamed);
				pOriginalProcessor[1]->getMesh()->addAttrMeshFeatures(*mfl2, StrAttrFeatureUnnamed);

				for (auto iter = vMatched.begin(); iter != vMatched.end(); ++iter)
				{
					if (iter->m_idx1 == iter->m_idx2) iter->m_note = 1;
					else iter->m_note = -1;
				}
				mpSwap = vMatched;
			}
			break;
		}

		/*-------------------------------------------------------------*/
// 		timer.startTimer();
// 		for (int group = 0; group < groupNum; ++group)
// 		{
// 			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
// 			double vPara[] = {40, 0.7};
// 			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 1, vPara, cout);
// 			ofs << "  Match results(method=" << 1 << ",t=" << 40 << "): "; 
// 			printVectorMatchPair(vMatched, ofs);
// 			validMatched = countValidMatchPair(vMatched);
// 			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
// 			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
// 			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
// 		}
// 		timer.stopTimer();
// 		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
// 		test_cat++;

		/*-------------------------------------------------------------*/
// 		timer.startTimer();
// 		for (int group = 0; group < groupNum; ++group)
// 		{
// 			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
// 			double vPara[] = {80, 0.7};
// 			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 1, vPara, cout);
// 			ofs << "  Match results(method=" << 1 << ",t=" << 80 << "): "; 
// 			printVectorMatchPair(vMatched, ofs);
// 			validMatched = countValidMatchPair(vMatched);
// 			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
// 			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
// 			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
// 		}
// 		timer.stopTimer();
// 		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
// 		test_cat++;

		/*-------------------------------------------------------------*/

// 		timer.startTimer();
// 		for (int group = 0; group < groupNum; ++group)
// 		{
// 			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
// 			double vPara[] = {40, 0.7};
// 			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 0, vPara, cout);
// 			ofs << "  Match results(method=" << 1 << ",t=" << 40 << "): "; 
// 			printVectorMatchPair(vMatched, ofs);
// 			validMatched = countValidMatchPair(vMatched);
// 			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
// 			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
// 			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
// 		}
// 		timer.stopTimer();
// 		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
// 
// 		test_cat++;

		/*-------------------------------------------------------------*/
// 		if (N <= 30)
// 		{
// 			timer.startTimer();
// 			for (int group = 0; group < groupNum; ++group)
// 			{
// 				const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
// 				double vPara[] = {0.52, 0.02, 0.01};
// 				PairMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 0, vPara, cout, false);
// 				ofs << "  Match results(method=PairGraph: "; 
// 				printVectorMatchPair(vMatched, ofs);
// 				validMatched = countValidMatchPair(vMatched);
// 				sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
// 				if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
// 				else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
// 			}
// 			timer.stopTimer();
// 			cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
// 		}	
// 		test_cat++;

		/*-------------------------------------------------------------*/
// 		timer.startTimer();
// 		for (int group = 0; group < groupNum; ++group)
// 		{
// 			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
// 			double vPara[] = {10, 1280, 0.52};
// 			// 			vector<double> vTimes;
// 			// 			for (int i = 0; i < 8; ++i) vTimes.push_back(10 * pow(2.,i));
// 			// 			SimplePointMatching(pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vTimes, vMatched, false);
// 			HKSMatchingExt(pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 0, vPara, cout, false);
// 			ofs << "  Match results(method=Simple: "; 
// 			printVectorMatchPair(vMatched, ofs);
// 			validMatched = countValidMatchPair(vMatched);
// 			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
// 			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
// 			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
// 		}
// 		timer.stopTimer();
// 		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
// 		test_cat++;
		/*-------------------------------------------------------------*/

		/*-------------------------------------------------------------*/
		timer.startTimer();
		for (int group = 0; group < groupNum; ++group)
		{
			const vector<int>& vinput1 = vvInput1[group], &vinput2 = vvInput2[group];
			int anchor = double(rand())/double(RAND_MAX) * vsize;
			double vPara[] = {10, 1280, 0.52, anchor};
			// 			vector<double> vTimes;
			// 			for (int i = 0; i < 8; ++i) vTimes.push_back(10 * pow(2.,i));
			// 			SimplePointMatching(pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vTimes, vMatched, false);
			HKSMatchingExt(pOriginalProcessor[0], pOriginalProcessor[1], vinput1, vinput2, vMatched, 1, vPara, cout, false);
			ofs << "  Match results(method=Simple: "; 
			printVectorMatchPair(vMatched, ofs);
			validMatched = countValidMatchPair(vMatched);
			sum_ratio_matched[test_cat] += (double)validMatched/(double)N;
			if (vMatched.empty()) sum_ratio_correct[test_cat] += 0.;
			else sum_ratio_correct[test_cat] += (double)validMatched/(double)vMatched.size();
		}
		timer.stopTimer();
		cout << "Average time: " << timer.getElapsedTime()/groupNum << endl;
		test_cat++;
		/*-------------------------------------------------------------*/

		ofs2 << int(N * (1.+outlier_ratio));
		for (int i = 0; i < test_cat; ++i)
		{
			ofs2 << ", " << sum_ratio_matched[i] / (double)groupNum << ", " << sum_ratio_correct[i] / (double)groupNum;
		}
		ofs2 << endl;
	}
}

void ShapeMatcher::localCorrespondenceTesting()
{
	assert(pOriginalMesh[0]->vertCount() == pOriginalMesh[1]->vertCount());
	int vsize = pOriginalMesh[0]->vertCount();

	ofstream ofs("output/evaluate_hkt.txt");
	ofstream ofs2("output/result_hkt.csv");

	const int groupNum = 30;
	const int anchorNum = 20;
	vector<int> vTest(groupNum);
	generate(vTest.begin(), vTest.end(), [&](){return double(rand())/double(RAND_MAX)*vsize;});

	vector<MatchPair> vAnchors;
	for (int i = 0; i < anchorNum; ++i)
	{
		int a = double(rand())/double(RAND_MAX)*vsize;
		vAnchors.push_back(MatchPair(a, a));
	}

	vector<MatchPair> matched;
	int valid_match;

	double match_ratio[20] = {0.};
	double correct_ratio[20] = {0.};
	int nCat;
	for (int n = 0; n < groupNum; ++n)
	{
		const int vi = vTest[n];
		vector<int> vNeighbor1 = pOriginalMesh[0]->getVertNeighborVerts(vi, 1, false);
		vector<int> vNeighbor2 = pOriginalMesh[1]->getVertNeighborVerts(vi, 1, false);
		if (vNeighbor1.size() != vNeighbor2.size())
		{
			cout << "Not compatible neighbors!";
		}
		else cout << "Neighbor of " << vi << ", size = " << vNeighbor1.size() << endl;
		nCat = 0;

		{
			double vPara[] = {5, 0.05};
			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, 0, vPara, cout, false);
			valid_match = countValidMatchPair(matched);
			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
			if (matched.empty()) correct_ratio[nCat] += 0;
			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
			++nCat;
		}
// 
// 
// 		{
// 			double vPara[] = {5, 0.5, vi, vi};
// 			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, 2, vPara, cout, false);
// 			valid_match = countValidMatchPair(matched);
// 			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
// 			if (matched.empty()) correct_ratio[nCat] += 0;
// 			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
// 			++nCat;
// 		}
// 
		{
			double vPara[] = {5, 0.6, 15};
			TensorMatchingExt(m_ep, pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, 3, vPara, cout, false);
			valid_match = countValidMatchPair(matched);
			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
			if (matched.empty()) correct_ratio[nCat] += 0;
			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
			++nCat;
		}

// 		{
// 			double vPara[] = {5, 10, 0.5};
// 			HKSMatchingExt(pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, 0, vPara, cout, false);
// 			valid_match = countValidMatchPair(matched);
// 			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
// 			if (matched.empty()) correct_ratio[nCat] += 0;
// 			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
// 			++nCat;
// 		}
		
		{
			double vPara[] = {5, 20, 0.5, vi, vi};
			HKSMatchingExt(pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, 1, vPara, cout, false);
			valid_match = countValidMatchPair(matched);
			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
			if (matched.empty()) correct_ratio[nCat] += 0;
			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
			++nCat;
		}

// 		{
// 			HKCMatching(pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, vector<MatchPair>(vAnchors.begin(), vAnchors.begin()+10), 40, 0.1);
// 			valid_match = countValidMatchPair(matched);
// 			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
// 			if (matched.empty()) correct_ratio[nCat] += 0;
// 			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
// 			++nCat;
// 		}

		{
			HKCMatching(pOriginalProcessor[0], pOriginalProcessor[1], vNeighbor1, vNeighbor2, matched, vAnchors, 40, 0.1);
			valid_match = countValidMatchPair(matched);
			match_ratio[nCat] += (double)valid_match / (double)vNeighbor1.size();
			if (matched.empty()) correct_ratio[nCat] += 0;
			else correct_ratio[nCat] += (double)valid_match / (double)matched.size();
			++nCat;
		}
	}

	for (int i = 0; i < nCat; ++i)
	{
		cout << match_ratio[i] / (double)groupNum << ", ";// << sum_ratio_correct[i] / (double)groupNum;
	}
	cout << endl;
	for (int i = 0; i < nCat; ++i)
	{
		cout << correct_ratio[i] / (double)groupNum << ", ";// << sum_ratio_correct[i] / (double)groupNum;
	}
	ofs2 << endl;
	
}

void ShapeMatcher::HKCMatching(MeshHelper* pmp1, MeshHelper* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, const std::vector<MatchPair>& vAnchorPair, double t, double thresh )
{
	vMatchedPair.clear();

	int vsize1 = vFeatures1.size(), vsize2 = vFeatures2.size();
	int anchorSize = vAnchorPair.size();
	vector<VecNd> vSig1(vsize1), vSig2(vsize2);

	{
		for (int i = 0; i < vsize1; ++i) 
		{
			vSig1[i].resize(anchorSize);
			for (int j = 0; j < anchorSize; ++j)
				vSig1[i][j] = calHK(pmp1, vAnchorPair[j].m_idx1, vFeatures1[i], t) /*/ vTrace1[j]*/;
		}
		for (int i = 0; i < vsize2; ++i) 
		{
			vSig2[i].resize(anchorSize);
			for (int j = 0; j < anchorSize; ++j)
				vSig2[i][j] = calHK(pmp2, vAnchorPair[j].m_idx2, vFeatures2[i], t) /*/ vTrace2[j]*/;
		}
	}

	double *hksSim = new double[vsize1 * vsize2];

	for (int i = 0; i < vsize1; ++i)
	{
		for (int j = 0; j < vsize2; ++j)
			hksSim[i * vsize2 + j] = vSig1[i].distEculidean2(vSig2[j]) / anchorSize;
	}

	double vMin = 1e10 - 1;
	while (vMin < 1e10 - 0.1)
	{
		double *pMin = min_element(hksSim, hksSim + vsize1 * vsize2);
		vMin = *pMin;

		if (vMin > thresh) break;

		size_t pos = pMin - hksSim;
		int i2 = pos % vsize2;
		int i1 = pos / vsize2;
		vMatchedPair.push_back(MatchPair(vFeatures1[i1], vFeatures2[i2], vMin));

		for (int k = 0; k < vsize1; ++k)
			hksSim[k*vsize2 + i2] = 1e10;
		for (int k = 0; k < vsize2; ++k)
			hksSim[i1*vsize2 + k] = 1e10;
	}

	delete []hksSim;
}

void ShapeMatcher::generateExampleMatching( int n )
{
	if (g_configMgr.getConfigValueInt("GROUND_TRUTH_AVAILABLE") != 1)
		return;

	int vsize = pOriginalMesh[0]->vertCount();

	const int anchorNum = n;
	vector<MatchPair> vAnchors;
	for (int i = 0; i < anchorNum; ++i)
	{
		int a = double(rand())/double(RAND_MAX)*vsize;
		vAnchors.push_back(MatchPair(a, a));
	}
	
	forceInitialAnchors(vAnchors);
}

void ShapeMatcher::loadGroundTruth( const std::string& filename )
{
	ifstream fin(filename.c_str());
	if (!fin) {
		std::cerr << "Ground truth file not available" << std::endl;
		return;
	}

	mMatchGroundTruth.clear();

	int size;
	fin >> size;
	for (int i = 0; i < size; ++i) {
		int v1, v2;
		fin >> v1 >> v2;
		mMatchGroundTruth.insert(std::make_pair(v1, v2));
	}

	m_bHasGroundTruth = true;
}

void ShapeMatcher::autoGroundTruth()
{
	assert(pOriginalMesh[0]->vertCount() == pOriginalMesh[1]->vertCount());
	int meshSize = pOriginalMesh[0]->vertCount();
	mMatchGroundTruth.clear();
	for (int i = 0; i < meshSize; ++i)
		mMatchGroundTruth.insert(std::make_pair(i, i));

	m_bHasGroundTruth = true;
}
