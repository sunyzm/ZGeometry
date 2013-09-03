#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <utility>
#include <engine.h>
#include "DifferentialMeshProcessor.h"
#include "ZMesh/MeshPyramid.h"

class MatchPair
{
public:
	int		m_idx1;
	int		m_idx2;
	double	m_tl;		//lower time
	double	m_tu;		//upper time
	int		m_tn;
	double  m_score;
	int		m_note;		

public:
	MatchPair() { m_idx1 = -1; m_idx2 = -1; }
	MatchPair(int i1, int i2, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_score = score; m_tl = 0.0; m_tn = 0;  m_note = 0;}
	MatchPair(int i1, int i2, double tl, int tn, double score = 0) { m_idx1 = i1; m_idx2 = i2; m_tl = tl; m_tn = tn; m_tu = m_tl * pow(2., tn-1); m_score = score; m_note = 0;}
	operator std::pair<int, int>() const { return std::make_pair(m_idx1, m_idx2); }
	
	friend bool operator== (const MatchPair& mp1, const MatchPair& mp2) { return (mp1.m_idx1 == mp2.m_idx1 && mp1.m_idx2 == mp2.m_idx2); }
 	friend bool operator< (const MatchPair& mp1, const MatchPair& mp2)
 	{
 		return mp1.m_idx1 < mp2.m_idx1 || (mp1.m_idx1 == mp2.m_idx1 && mp1.m_idx2 < mp2.m_idx2);
 	}
 
 	friend bool operator> (const MatchPair& mp1, const MatchPair& mp2)
 	{
 		return mp2 < mp1;
 	}

	static std::vector<std::pair<int, int> > ToPairVector(const std::vector<MatchPair>& vmp)
    {
		std::vector<std::pair<int, int> > vp;
		for (auto iter = vmp.begin(); iter != vmp.end(); ++iter) 
			vp.push_back( (std::pair<int, int>)(*iter) );
		return vp;
	}
};

class PairScoreCompare
{
public:
	bool operator()(const MatchPair& Left, const MatchPair& Right) const { return ( Left.m_score > Right.m_score ); }
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

class HKSFeature : public MeshFeature
{
public:
	double	m_tl;	    // lower time
	double	m_tu;	    // upper time
	int		m_tn;	    // num of scales; log_2(tu/tl)
	int     minOrMax;	// -1:min; 1:max
public:
	HKSFeature() : MeshFeature() { m_index = 0; m_scale = -1; m_tu = m_tl = 0.0; m_tn = 0; }
	HKSFeature(int i, int s, int state = 1) : minOrMax(state) { m_index = i; m_scale = s; m_tu = m_tl = 0.0; m_tn = 0; }
	HKSFeature(int x, double tl, double tu) { m_index = x; m_tl = tl; m_tu = tu; m_tn = 0; }
	HKSFeature(int x, double tl, int n) { m_index = x; m_tn = n; m_tu = tl; }
	HKSFeature(int x, double tl, double tu, int n) { m_index = x; m_tn = n; m_tl = tl; m_tu = tu; }
	HKSFeature(int x, double tl, double tu, int n, int scale) { m_index = x; m_tn = n; m_tl = tl; m_tu = tu; m_scale = scale; }
	HKSFeature(const HKSFeature& f) { m_index = f.m_index; m_scale = f.m_scale; m_tl = f.m_tl; m_tu = f.m_tu; m_tn = f.m_tn; minOrMax = f.minOrMax; }
	void setTimes(double tl, double tu, double tn) { m_tl = tl; m_tu = tu; m_tn = tn; }
};

class PointParam : public VectorND
{
public:
	PointParam& operator =(const PointParam& hkp);
	double m_votes;
	virtual void clear();
};

class ParamManager
{
public:
	enum ParamType {HKParam, BHParam};
	int para_dim;
	std::vector<PointParam> vCoord;
	std::vector<PointParam> vSignature;
	void initialize(DifferentialMeshProcessor* p) { pMP = p; }
	const std::vector<PointParam>& getParam() const { return vCoord; }

	DifferentialMeshProcessor* pMP;

	void computeHKParam(const std::vector<int>& anchors, double t = 30.0);
	void computeBHParam(const std::vector<int>& anchors);
	void para_computeHKC(const std::vector<int>& anchors, double t = 30.0);
	void para_computeHKS(const std::vector<double>& times);
};

class DiffusionShapeMatcher
{
public:
	std::vector<MatchPair> mpSwap;
	void swapMP() { std::swap(m_vFeatureMatchingResults[m_nAlreadyMatchedLevel], mpSwap); }
public:
	DiffusionShapeMatcher();
	~DiffusionShapeMatcher();

	/* core functions */
	void    initialize(DifferentialMeshProcessor* pMP1, DifferentialMeshProcessor* pMP2, Engine *ep);
	void	constructPyramid(int n, double ratio, std::ostream& ostr);
	void	detectFeatures(int obj, int ring = 2, int scale = 1, double tvalue = DEFAULT_FEATURE_TIMESCALE, double talpha = DEFAULT_T_MULTIPLIER, double thresh = DEFAULT_EXTREAMA_THRESH);
	void	matchFeatures(std::ostream& flog, double matchThresh = DEFAULT_MATCH_THRESH);
	void    matchFeatureSimple();
	void	matchFeaturesTensor_deprecate(std::ostream& flog, double timescale, double thresh);
	void	refineRegister(std::ostream& flog);
	void    refineRegister2(std::ostream& flog);
	
	void	evaluateRegistration();

	/* testing functions */
	void    registerTesting1();
	void    regsiterTesting2();
	void    dataTesting1();
	void    sparseMatchingTesting();
	void    generateExampleMatching(int n); 
	void    localCorrespondenceTesting();

	/* attributes access */
	DifferentialMeshProcessor* getMeshProcessor(int obj, int level) { return liteMP[obj].at(level); }
	int     getPyramidLevels() const { return m_nPyramidLevels; }
	bool	isPyramidBuilt() const { return m_bPyramidBuilt; }
	int		getTotalRegistrationLevels() const { return m_nRegistrationLevels; }
	void	setRegistrationLevels(int val);
	void	setEngine(Engine* ep) { m_ep = ep; }
	const MeshPyramid& getMeshPyramid(int obj) const { return meshPyramids[obj]; }
	CMesh*  getMesh(int obj, int level = 0) const;
	int     getAlreadyRegisteredLevel() const { return m_nAlreadyRegisteredLevel; }
	int		getAlreadyMatchedLevel() const { return m_nAlreadyMatchedLevel; }
	const std::vector<MatchPair>& getInitialMatchedFeaturePairs() const;
	const std::vector<MatchPair>& getMatchedFeaturesResults(int level) const;
	const std::vector<MatchPair>& getRegistrationResults(int level) const;
	std::vector<HKSFeature>& getSparseFeatures(int obj) { return m_vFeatures[obj]; }
	const std::vector<HKSFeature>& getSparseFeatures_const(int obj) const { return m_vFeatures[obj]; }
	bool	loadInitialFeaturePairs(const std::string& filename);
	void	forceInitialAnchors(const std::vector<MatchPair>& mp);
	void	loadGroundTruth(const std::string& filename);
	bool	hasGroundTruth() const { return m_bHasGroundTruth; }

	int		id2Index(int obj, int vid, int level) const { return meshPyramids[obj].m_Id2IndexMap[vid][level]; }
	void    dumpIndexMap(const std::string& filename) const;
	void	readInRandPair(const std::string& filename);
	
	void evaluateWithGroundTruth(const std::vector<MatchPair>& vIdMatchPair);
	static double evaluateDistortion(const std::vector<MatchPair>& vIdMatchPair, const CMesh* mesh1, const CMesh* mesh2, const std::vector<std::pair<double, double> >& vRandPair, int rand_start = 0);
	static double evaluateDistance(const DifferentialMeshProcessor& mp1, const DifferentialMeshProcessor& mp2, DistanceType distType, const std::vector<double>& vParam, const std::vector<std::pair<double, double> >& vRandPair, int rand_start = 0);

	std::vector<std::pair<double, double> > m_randPairs;

	// static constants
	static const double DEFAULT_C_RATIO;
	static const double DEFAULT_RANK_EPSILON;
	static const double SPARSIFY_EPSILON;
	static const double DEFAULT_FEATURE_TIMESCALE;
	static const double DEFAULT_T_MULTIPLIER;
	static const double DEFAULT_MATCH_TIME_LOW;
	static const double DEFAULT_MATCH_THRESH;
	static const double DEFAULT_EXTREAMA_THRESH;
	static const double DEFAULT_REGISTER_TIMESCALE;
	static const int	DEFAULT_PYRAMID_LEVELS;
	static const int	MAXIMAL_PYRAMID_LEVELS;
	static const int    NUM_OF_EIGVAL_FOR_ESTIMATE;

private:
	Engine* m_ep;
	CMesh* pOriginalMesh[2];
	DifferentialMeshProcessor* pOriginalProcessor[2];
	MeshPyramid meshPyramids[2];
	std::vector<DifferentialMeshProcessor*> liteMP[2];
	std::vector<std::vector<HKSFeature> > m_vFeatures;	// original detected fine features

	std::vector<std::vector<MatchPair> > m_vFeatureMatchingResults;
	std::vector<std::vector<MatchPair> > m_vRegistrationResutls;

	bool			   m_bHasGroundTruth;
	std::map<int, int> mMatchGroundTruth;

	bool					m_bInitialized;
	bool					m_bPyramidBuilt;
	bool					m_bFeatureDetected;
	bool					m_bFeatureMatched;
	int						m_nPyramidLevels;           // >= 1
	int						m_nRegistrationLevels;      // [1,..,m_nPyramidLevels]
	int						m_nAlreadyRegisteredLevel;	// [0,..,m_nRegistrationLevels]
	int						m_nAlreadyMatchedLevel;		// [1,...,m_nRegistrationLevels]
	int						m_nBaseEigensMatch, m_nBaseEigensRegister;
	double					m_registerTimescale;
	ParamManager			m_ParamMgr[2];

public:
	/* helper functions */
	static void	calVertexSignature( const DifferentialMeshProcessor* pOriginalProcessor, const HKSFeature& hf, VectorND& sig );
	static void SimplePointMatching(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, const std::vector<double>& vTimes, std::vector<MatchPair>& machedResult, bool verbose = false);
	static void PairGraphMatching(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<HKSFeature>& vFeatures1, const std::vector<HKSFeature>& vFeatures2, std::vector<MatchPair>& vMatchedPair, double para_thresh, bool verbose = false);
	
	static void HKSMatchingExt(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, int method, double vPara[], std::ostream& logout, bool verbose = false);
	static double PairMatchingExt(Engine* ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, int method, double vPara[], std::ostream& logout, bool verbose = false);
	static double TensorMatchingExt(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPairs, int highOrderFeatureType, double vPara[], std::ostream& logout, bool verbose = false);
	static void HKCMatching(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& vMatchedPair, const std::vector<MatchPair>& vAnchorPair, double t, double thresh);

	static double TensorGraphMatching6(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& matched, double para_t, double para_thresh, bool verbose = false);
	static double TensorGraphMatching6(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<HKSFeature>& vFeatures1, const std::vector<HKSFeature>& vFeatures2, std::vector<MatchPair>& matched, double para_t, double para_thresh, bool verbose = false);
	static double TensorGraphMatching3(Engine *ep, const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, const std::vector<int>& vFeatures1, const std::vector<int>& vFeatures2, std::vector<MatchPair>& matched, double t, double thresh);
	static void ComputeTensorFeature12( const DifferentialMeshProcessor* pmp, int i, int j, int k, double t1, double t2, double* sang );
	static void ComputeTensorFeature6( const DifferentialMeshProcessor* pmp, int i, int j, int k, double t, double* sang, bool sweep = false);
	static void ComputeTensorFeature3( const DifferentialMeshProcessor* pmp, int i, int j, int k, double t, double* sang);
	static void ComputeTensorFeatureAnchor(const DifferentialMeshProcessor* pmp, int i, int j, int k, int origin, double t, double* sang);
	void    prepareHeatRegistration( double regTime );
	double computeMatchScore(int idx1, int idx2, double sigma = 0.02) const;
	int		searchVertexMatch( const int vt, const int vj, const int level, const int ring, double& score, int uppper_level = -1 );
	void    getVertexCover(int obj, int vidx, int level, int upper_level, int ring,  std::vector<int>& vCoveredIdx) const;
	static double  calPointHksDissimilarity(const DifferentialMeshProcessor* pmp1, const DifferentialMeshProcessor* pmp2, int i1, int i2, const std::vector<double>& vTimes, int mode = 0);
	void autoGroundTruth();
};