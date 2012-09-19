#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_
#include <complex>
#include <vector>
#include <string>
#include <engine.h>
#include "Mesh.h"

#define LBO_EPS 0
#define LBO_COT 1
#define LBO_GRAPH 2

// const double PZERO  = 0.00001;
// const double NZERO  = -0.00001;
// const int HOLE_SIZE = 10;

class ManifoldHarmonics;

class Laplacian
{
public:
	int size;
	std::vector<double> vWeights;
	std::vector<int> vII, vJJ;
	std::vector<double> vSS;

	enum LaplacianType {Graph, CotFormula};
	void computeLaplacian(const CMesh* tmesh, LaplacianType laplacianType = CotFormula);
	void decompose(ManifoldHarmonics& mhb, int nEig, Engine *ep) const;
	double innerProduct(const std::vector<double>& vf, const std::vector<double>& vg) const;
};

class ManifoldBasis
{
public:
	std::vector<double> m_vec;	// eigenvector
	double m_val;			    // eigenvalue
};

class ManifoldHarmonics
{
public:
	std::vector<ManifoldBasis> m_func;	// manifold harmonic basis
	int m_size;	    // shape size
	int m_nEigFunc; // number of basis < shape size
public:
	ManifoldHarmonics() : m_size(0), m_nEigFunc(0) {}
	bool decompLaplacian(Engine *ep, const CMesh *tmesh, int nEigFunc, short lbo_type = LBO_COT);
	void write(const std::string& meshPath) const;
	void read(const std::string& meshPath);
};

#endif