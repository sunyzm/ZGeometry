#pragma once
#include "DifferentialMeshProcessor.h"
#include "MeshPyramid.h"
#include <vector>

struct MatchPair
{
	MatchPair(int i1, int i2) : m_idx1(i1), m_idx2(i2), m_matchScore(0) {}
	MatchPair(int i1, int i2, double s) : m_idx1(i1), m_idx2(i2), m_matchScore(s) {}

	int m_idx1, m_idx2;
	double m_matchScore;
};

class DiffusionShapeMapper
{
public:
	DiffusionShapeMapper() { pMP[0] = pMP[1] = NULL; pOriginalMesh[0] = pOriginalMesh[1] = NULL; }
	void initialize(const DifferentialMeshProcessor* pMP1, const DifferentialMeshProcessor* pMP2);
	void buildPyramid();
	void detectFeatures();
	void matchFeatures();
	const std::vector<MatchPair>& getMatchedFeatures();

private:
	CMesh* pOriginalMesh[2];
	DifferentialMeshProcessor *pMP[2];
	MeshPyramid mPyramid[2];

	std::vector<int> vFeatures1, vFeatures2;

};