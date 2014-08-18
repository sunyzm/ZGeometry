#include "MeshLaplacian.h"
#include <ctime>
#include <algorithm>
#include <cstdio>
#include <tuple>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <ZGeom/arithmetic.h>
#include <ZGeom/EigenCompute.h>
#include <ZGeom/MeshProcessing.h>
#include <ZGeom/zassert.h>

using namespace std;
using ZGeom::PI;
using ZGeom::uint;

void MeshLaplacian::constructTutte( const CMesh* tmesh )
{
	mOrder = tmesh->vertCount();
	std::vector<std::tuple<int,int,double> > vElem;
	std::vector<double> vWeights(mOrder);

	for (int i = 0; i < mOrder; ++i) {
		const CVertex* vi = tmesh->getVertex(i);
		vector<int> vNeighbors;
		tmesh->vertRingNeighborVerts(i, 1, vNeighbors, false);
		const int valence = vNeighbors.size();

		for (int j = 0; j < valence; ++j) {
			vElem.push_back(std::make_tuple(i+1, vNeighbors[j]+1, 1.0));
		}
		vElem.push_back(std::make_tuple(i+1, i+1, -double(valence)));	// diagonal elements

		vWeights[i] = double(valence);
	}

	mLS.convertFromCOO(mOrder, mOrder, vElem);
	mW.convertFromDiagonal(vWeights);

	mLaplacianType = Tutte;
	mSymmetric = false;
}

void MeshLaplacian::constructUmbrella( const CMesh* tmesh )
{
	mOrder = tmesh->vertCount();
	ConstructMeshMatrix(*tmesh, ZGeom::MM_GRAPH_LAPLACE, mLS);
	mLS.scale(-1);
	mW.setToIdentity(mOrder);	// set vertex weight matrix to identity; the attained Laplacian becomes symmetric

	mLaplacianType = Umbrella;
	mSymmetric = true;
}

void MeshLaplacian::constructNormalizedUmbrella( const CMesh* tmesh )
{
	/* L = D^(-1/2) * (D-A) * D^(-1/2) */
	mOrder = tmesh->vertCount();
	ConstructMeshMatrix(*tmesh, ZGeom::MM_NORMALIZED_GRAPH_LAPLACE, mLS);
	mLS.scale(-1);
	mW.setToIdentity(mOrder);

	mLaplacianType = NormalizedUmbrella;
	mSymmetric = true;
}

/* Construct negative discrete Laplace operator */
void MeshLaplacian::constructCotFormula( const CMesh* tmesh )
{
	const int vertCount = tmesh->vertCount();
	mOrder = vertCount;
	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights(vertCount);

	std::vector<double> diagW(vertCount, 0);
	for (int vIdx = 0; vIdx < vertCount; ++vIdx)  //for every vertex
	{
		double amix = 0; //mixed area
		const CVertex* vi = tmesh->getVertex(vIdx);
		for (int j = 0; j < vi->outValence(); ++j) {
			const CHalfEdge* e0 = vi->getHalfEdge(j);
			const CHalfEdge* e1 = e0->nextHalfEdge();
			const CHalfEdge* e2 = e1->nextHalfEdge();
			const CVertex* vj = e0->vert(1);
			const CHalfEdge* e2twin = e2->twinHalfEdge();

			double len0 = e0->getLength();
			double len1 = e1->getLength();
			double len2 = e2->getLength();
			amix += ZGeom::calMixedTriArea(len0, len1, len2);
			
			double cot_a11(0), cot_a12(0), cot_a21(0), cot_a22(0);
			ZGeom::triangleCot(len0, len1, len2, cot_a11, cot_a12);

			if (e0->twinHalfEdge() != NULL) {
				const CHalfEdge* e10 = e0->twinHalfEdge();
				const CHalfEdge* e11 = e10->nextHalfEdge();
				const CHalfEdge* e12 = e11->nextHalfEdge();
				len1 = e11->getLength();
				len2 = e12->getLength();

				ZGeom::triangleCot(len0, len1, len2, cot_a21, cot_a22);
			}
			double cota = (cot_a11 + cot_a21) / 2.0;

			vII.push_back(vIdx + 1);
			vJJ.push_back(vj->getIndex() + 1);
			vSS.push_back(cota);
			diagW[vIdx] -= cota;

			if (e2twin == NULL) { //met an boundary fan edge
				const CHalfEdge* e20 = e2;
				const CHalfEdge* e21 = e20->nextHalfEdge();
				const CHalfEdge* e22 = e21->nextHalfEdge();
				len0 = e20->getLength();
				len1 = e21->getLength();
				len2 = e22->getLength();
				double cot_a1, cot_a2;
				ZGeom::triangleCot(len0, len1, len2, cot_a1, cot_a2);
				double cota = cot_a1 / 2.0;
				
				vII.push_back(vIdx + 1);
				vJJ.push_back(e20->getVertIndex(0) + 1);
				vSS.push_back(cota);
				diagW[vIdx] -= cota;
			}
		} // for each incident halfedge

		vWeights[vIdx] = amix;
	}

	for (int vIdx = 0; vIdx < vertCount; ++vIdx) {
		vII.push_back(vIdx + 1);
		vJJ.push_back(vIdx + 1);
		vSS.push_back(diagW[vIdx]);
	}

	mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);

	mLaplacianType = CotFormula;
	mSymmetric = false;
}

void MeshLaplacian::constructSymCot( const CMesh* tmesh )
{
	const int vertCount = tmesh->vertCount();
	constructCotFormula(tmesh);	
	mW.setToIdentity(vertCount);

	mLaplacianType = SymCot;
	mSymmetric = true;
}

void MeshLaplacian::constructAnisotropic1( const CMesh* tmesh )
{
	using namespace std;
	const int vertCount = mOrder = tmesh->vertCount();
	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals();
	const std::vector<double>& vMeanCurvatures = tmesh->getMeanCurvature();
	const std::vector<double>& vMixedAreas = tmesh->getVertMixedAreas();
	vector<std::tuple<int,int,double> > vSparseElements;

	int nRing = 1;
	double hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	double hPara2 = 0.2;

	std::vector<double> vDiag(vertCount, 0);
	
	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->getVertex(vi); 
		int outValence = pvi->outValence();
		for (int j = 0; j < outValence; ++j) {
			const CHalfEdge* phe = pvi->getHalfEdge(j);
		    const CVertex* pvj = phe->vert(1);
			int vj = pvj->getIndex();
			if (vi > vj && phe->twinHalfEdge() != NULL) continue;

			Vector3D vij = pvj->getPosition() - pvi->getPosition();
			double w1 = std::exp(-vij.length2() / hPara1);
			double w2 = std::exp(-pow(
				(fabs(dotProduct3D(vVertNormals[vi], vij)) + std::fabs(dotProduct3D(vVertNormals[vj], vij)))
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
	vWeights = vMixedAreas;
	mLS.convertFromCOO(vertCount, vertCount, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);
#ifdef _DEBUG
	std::cout << "Anisotropic Laplacian is symmetric? " << (bool)mLS.testSymmetric() << '\n';
#endif
	
	mLaplacianType = Anisotropic1;
	mSymmetric = false;
}

void MeshLaplacian::constructAnisotropic2(const CMesh* tmesh)
{
	using namespace std;
	const int vertCount = mOrder = tmesh->vertCount();
	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals();
	const std::vector<double>& vMeanCurvatures = tmesh->getMeanCurvature();
	const std::vector<double>& vMixedAreas = tmesh->getVertMixedAreas();
	vector<std::tuple<int,int,double> > vSparseElements;	
	double avgMeanCurv(0);
	for (double a : vMeanCurvatures) avgMeanCurv += fabs(a);
	avgMeanCurv /= vertCount;
	int nRing = 1;
	double hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	double hPara2 = 0.2;

	std::vector<double> vDiag(vertCount, 0);
	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->getVertex(vi); 
		int outValence = pvi->outValence();
		for (int j = 0; j < outValence; ++j) {
			const CHalfEdge* phe = pvi->getHalfEdge(j);
			const CVertex* pvj = phe->vert(1);
			int vj = pvj->getIndex();
			if (vi > vj && phe->twinHalfEdge() != NULL) continue;

			Vector3D vij = pvj->getPosition() - pvi->getPosition();
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
	vWeights = vMixedAreas;

	mLS.convertFromCOO(vertCount, vertCount, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);	

	mLaplacianType = Anisotropic2;
	mSymmetric = false;
}

void MeshLaplacian::constructAnisotropic3( const CMesh* tmesh, int nRing, double hPara1, double hPara2 )
{
	// based on curvature difference or bilateral filtering
	const int vertCount = mOrder = tmesh->vertCount();
	nRing = 1;
	hPara1 = 2 * std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);

	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals();
	const std::vector<double>& vMeanCurvatures = tmesh->getMeanCurvature();

	vector<std::tuple<int,int,double> > vSparseElements;

	for (int vi = 0; vi < vertCount; ++vi) {
		const CVertex* pvi = tmesh->getVertex(vi); 
		vector<int> vAdjFaces = tmesh->getVertexAdjacentFaceIdx(vi, nRing);
		for (int fi : vAdjFaces) {
			const CFace* pfi = tmesh->getFace(fi);
			double face_area = pfi->calArea();
			for (int k = 0; k < 3; ++k) {
				int vki = pfi->getVertexIndex(k);
				if (vi == vki) continue;
				const CVertex* pvk = pfi->getVertex(k);

				Vector3D vpw = pvi->getPosition() - pvk->getPosition();
//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				double w1 = std::exp(-vpw.length2() / hPara1);
//				double w2 = std::exp(-std::pow(vMeanCurvatures[vi] - vMeanCurvatures[vki], 2)  / hPara2);
				double w2 = std::exp(-std::pow(dotProduct3D(vVertNormals[vi], vpw), 2) / hPara2);

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

	mLaplacianType = Anisotropic1;
	mSymmetric = true;
}

void MeshLaplacian::constructAnisotropic4(const CMesh* tmesh, int ringT, double hPara1, double hPara2)
{
	// similar to bilateral filtering
	mOrder = tmesh->vertCount();

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;

	vector< std::tuple<int,int,double> > vSparseElements;

	ringT = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
//	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = tmesh->getAvgEdgeLength();
	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals();

	for (int vi = 0; vi < mOrder; ++vi)	{
		const CVertex* pvi = tmesh->getVertex(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFaceIdx(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi) {
			const CFace* pfi = tmesh->getFace(vFaces[fi]);
			double face_area = pfi->calArea();
			for (int k = 0; k < 3; ++k)	{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex(k);

				double w1 = 1., w2 = 1.;
//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				w1 = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / hPara1);
//				w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
//				w2 = std::exp(-std::pow(dotProduct3D(pvi->getNormal(), pvi->getPosition() - pvk->getPosition()), 2) / hPara2);
				w2 = std::exp(dotProduct3D(vVertNormals[vi], pvk->getPosition() - pvi->getPosition()) / hPara2);
				
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

void MeshLaplacian::constructFromMesh5( const CMesh* tmesh )
{
	mOrder = tmesh->vertCount();

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;
	vector<std::tuple<int,int,double> > vSparseElements;

	int ringT = 5;
	double hPara1 = std::pow(tmesh->getAvgEdgeLength()*2, 2);
	
	for (int vi = 0; vi < mOrder; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFaceIdx(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace(vFaces[fi]);
			double face_area = pfi->calArea();
			for (int k = 0; k < 3; ++k)
			{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex(k);

				double svalue = std::exp(-(pvi->getPosition() - pvk->getPosition()).length2() / (4 * hPara1));
				svalue = svalue * face_area / (3 * 4 * PI * hPara1 * hPara1);

				vSparseElements.push_back(make_tuple(vi, vki, svalue));
			}
		}
	}

	vector<double> vDiag(mOrder, 0.);
	for (auto iter = begin(vSparseElements); iter != end(vSparseElements); ++iter)
	{
		int ii, jj; double ss;
		std::tie(ii, jj, ss) = *iter;

		vII.push_back(ii+1);
		vJJ.push_back(jj+1);
		vSS.push_back(ss);

		vDiag[ii] += -ss;
	}

	for (int i = 0; i < mOrder; ++i)
	{
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(vDiag[i]);
	}

	for (int k = 0; k < vII.size(); ++k)
	{
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
	if (mLaplacianType == Tutte || mLaplacianType == CotFormula ||
		mLaplacianType == Anisotropic1 || mLaplacianType == Anisotropic2) {
		this->decompose(nActualEigen, eng, es, true);
	}
	else if (mLaplacianType == Umbrella || mLaplacianType == SymCot) {
		this->decompose(nActualEigen, eng, es, false);
	}
	else throw std::logic_error("The symmetry property of Laplacian cannot be determined!");
}
