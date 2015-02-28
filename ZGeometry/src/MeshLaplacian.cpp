#include "MeshLaplacian.h"
#include <ctime>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <ZGeom/ZGeom.h>
#include <ZGeom/EigenCompute.h>
#include <ZGeom/zassert.h>

using namespace std;
using ZGeom::PI;
using ZGeom::uint;

void MeshLaplacian::constructAnisotropic1(CMesh* tmesh)
{
	using namespace std;
	const int vertCount = mOrder = tmesh->vertCount();
	const std::vector<ZGeom::Vec3d>& vVertNormals = tmesh->getVertNormals();
    const std::vector<double>& vMeanCurvatures = ZGeom::getMeshMeanCurvatures(*tmesh);
    const std::vector<double> vMixedAreas = computeMeshVertArea(*tmesh, ZGeom::VA_MIXED_VORONOI);
	vector<std::tuple<int,int,double> > vSparseElements;

	int nRing = 1;
	double hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	double hPara2 = 0.2;

	std::vector<double> vDiag(vertCount, 0);
	
	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->vert(vi); 
		int outValence = pvi->outValence();
		for (int j = 0; j < outValence; ++j) {
			const CHalfEdge* phe = pvi->getHalfEdge(j);
		    const CVertex* pvj = phe->vert(1);
			int vj = pvj->getIndex();
			if (vi > vj && phe->twinHalfEdge() != NULL) continue;

			ZGeom::Vec3d vij = pvj->pos() - pvi->pos();
			double w1 = std::exp(-vij.length2() / hPara1);
			double w2 = std::exp(-pow(
				(fabs(dot(vVertNormals[vi], vij)) + std::fabs(dot(vVertNormals[vj], vij)))
				/vij.length(), 2) / hPara2);
			double wij = w1 * w2;

			vSparseElements.push_back(make_tuple(vi, vj, wij));
			vSparseElements.push_back(make_tuple(vj, vi, wij));
			vDiag[vi] -= wij;
			vDiag[vj] -= wij;
		}
	}
	
	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	for (auto elem : vSparseElements) {
		int ii = std::get<0>(elem);
		int jj = std::get<1>(elem);
		double ss = std::get<2>(elem);
		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);
	}
	for (int i = 0; i < vertCount; ++i) {
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	std::vector<double> vWeights(vertCount, 1.0);
	//vWeights = vMixedAreas;
	mLS.convertFromCOO(vertCount, vertCount, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);

#ifdef _DEBUG
	std::cout << "Anisotropic Laplacian is symmetric? " << (bool)mLS.testSymmetric() << '\n';
#endif
	
	mSymmetric = true;
}

void MeshLaplacian::constructAnisotropic2(CMesh* tmesh)
{
	using namespace std;
	const int vertCount = mOrder = tmesh->vertCount();
	const std::vector<ZGeom::Vec3d>& vVertNormals = tmesh->getVertNormals();
    const std::vector<double>& vMeanCurvatures = ZGeom::getMeshMeanCurvatures(*tmesh);
    const std::vector<double> vMixedAreas = computeMeshVertArea(*tmesh, ZGeom::VA_MIXED_VORONOI);
	vector<std::tuple<int,int,double> > vSparseElements;	
	double avgMeanCurv(0);
	for (double a : vMeanCurvatures) avgMeanCurv += fabs(a);
	avgMeanCurv /= vertCount;
	int nRing = 1;
	double hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	double hPara2 = 0.2;

	std::vector<double> vDiag(vertCount, 0);
	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->vert(vi); 
		int outValence = pvi->outValence();
		for (int j = 0; j < outValence; ++j) {
			const CHalfEdge* phe = pvi->getHalfEdge(j);
			const CVertex* pvj = phe->vert(1);
			int vj = pvj->getIndex();
			if (vi > vj && phe->twinHalfEdge() != NULL) continue;

			ZGeom::Vec3d vij = pvj->pos() - pvi->pos();
			double w1 = std::exp(-vij.length2() / hPara1);
			double curvDiff = vMeanCurvatures[vi] - vMeanCurvatures[vj];
			double w2 = std::exp(-curvDiff*curvDiff / hPara2);
			double wij = w1 * w2;

			vSparseElements.push_back(make_tuple(vi, vj, wij));
			vSparseElements.push_back(make_tuple(vj, vi, wij));
			vDiag[vi] -= wij;
			vDiag[vj] -= wij;
		}
	}

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	for (auto elem : vSparseElements) {
		int ii = std::get<0>(elem);
		int jj = std::get<1>(elem);
		double ss = std::get<2>(elem);
		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);
	}
	for (int i = 0; i < vertCount; ++i) {
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	std::vector<double> vWeights(vertCount, 1.0);
	//vWeights = vMixedAreas;

	mLS.convertFromCOO(vertCount, vertCount, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);	

	mSymmetric = true;
}

void MeshLaplacian::constructAnisotropic3(CMesh* tmesh, int nRing, double hPara1, double hPara2)
{
	// based on curvature difference or bilateral filtering
	const int vertCount = mOrder = tmesh->vertCount();
	nRing = 1;
	hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);

	const std::vector<ZGeom::Vec3d>& vVertNormals = tmesh->getVertNormals();
    const std::vector<double>& vMeanCurvatures = ZGeom::getMeshMeanCurvatures(*tmesh);

	vector<std::tuple<int,int,double> > vSparseElements;

	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->vert(vi); 
		vector<int> vAdjFaces = tmesh->getVertexAdjacentFaceIdx(vi, nRing);
		for (int fi : vAdjFaces) {
			const CFace* pfi = tmesh->getFace(fi);
			double face_area = pfi->calArea();
			for (int k = 0; k < 3; ++k) {
				int vki = pfi->vertIdx(k);
				if (vi == vki) continue;
				const CVertex* pvk = pfi->vert(k);

				ZGeom::Vec3d vpw = pvi->pos() - pvk->pos();
//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				double w1 = std::exp(-vpw.length2() / hPara1);
//				double w2 = std::exp(-std::pow(vMeanCurvatures[vi] - vMeanCurvatures[vki], 2)  / hPara2);
				double w2 = std::exp(-std::pow(dot(vVertNormals[vi], vpw), 2) / hPara2);

				double svalue = face_area * w1 * w2;

				bool alreadyExisted = false;
				for (auto iterElem = vSparseElements.rbegin(); iterElem != vSparseElements.rend(); ++iterElem) {
					if (std::get<0>(*iterElem) == vi && std::get<1>(*iterElem) == vki) {
						std::get<2>(*iterElem) += svalue;
						alreadyExisted = true;
						break;
					}
				}
				if (!alreadyExisted)
					vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;
	vector<double> vDiag(vertCount, 0.);

	for (auto elem : vSparseElements) {
		int ii = std::get<0>(elem);
		int jj = std::get<1>(elem);
		double ss = std::get<2>(elem);
		
		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < vertCount; ++i) {
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}


	vWeights.resize(mOrder, 1.0);	
	mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
	std::cout << "Anisotropic Laplacian is symmetric? " << mLS.testSymmetric() << '\n';
	mW.convertFromDiagonal(vWeights);

	mSymmetric = true;
}

void MeshLaplacian::constructAnisotropic4(CMesh* tmesh, int nRing, double hPara1, double hPara2)
{
	// similar to bilateral filtering
	mOrder = tmesh->vertCount();

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;

	vector< std::tuple<int,int,double> > vSparseElements;

	nRing = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
//	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = tmesh->getAvgEdgeLength();
	const std::vector<ZGeom::Vec3d>& vVertNormals = tmesh->getVertNormals();

	for (int vi = 0; vi < mOrder; ++vi)	{
		const CVertex* pvi = tmesh->vert(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFaceIdx(vi, nRing);
		for (int fi = 0; fi < vFaces.size(); ++fi) {
			const CFace* pfi = tmesh->getFace(vFaces[fi]);
			double face_area = pfi->calArea();
			for (int k = 0; k < 3; ++k)	{
				int vki = pfi->vertIdx(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->vert(k);

				double w1 = 1., w2 = 1.;
//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				w1 = std::exp(-(pvi->pos() - pvk->pos()).length2() / hPara1);
//				w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
//				w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
				w2 = std::exp(dot(vVertNormals[vi], pvk->pos() - pvi->pos()) / hPara2);
				
				double svalue = w1 * w2;
				svalue *= face_area;

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(mOrder, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter) {
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < mOrder; ++i) {
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k) {
		vSS[k] /= vDiag[vII[k]-1];
	}

	vWeights.resize(mOrder, 1.0);	
	mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);

}

void MeshLaplacian::meshEigenDecompose(int nEig, ZGeom::MatlabEngineWrapper* eng, ZGeom::EigenSystem& es) const
{
	int nActualEigen = nEig;
	if (nEig == -1 || nEig >= mOrder) nActualEigen = mOrder - 1;
	if (mSymmetric) this->decompose(nActualEigen, eng, es, false);
	else this->decompose(nActualEigen, eng, es, true);
}

void MeshLaplacian::constructUmbrella(CMesh* tmesh)
{
    Laplacian::constructUmbrella(tmesh);
    mLaplacianType = Umbrella;
}

void MeshLaplacian::constructGeometricUmbrella(CMesh *tmesh)
{
    Laplacian::constructGeometricUmbrella(tmesh);
    mLaplacianType = Geometric;
}

void MeshLaplacian::constructNormalizedUmbrella(CMesh* tmesh)
{
    Laplacian::constructNormalizedUmbrella(tmesh);
    mLaplacianType = NormalizedUmbrella;
}

void MeshLaplacian::constructTutte(CMesh* tmesh)
{
    Laplacian::constructTutte(tmesh);
    mLaplacianType = Tutte;   
}

void MeshLaplacian::constructCotFormula(CMesh* tmesh)
{
    Laplacian::constructCotFormula(tmesh);
    mLaplacianType = CotFormula;
}

void MeshLaplacian::constructSymCot(CMesh* tmesh)
{
    Laplacian::constructSymCot(tmesh);
    mLaplacianType = SymCot;
}
