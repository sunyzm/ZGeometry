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

using namespace std;
using ZGeom::PI;
using ZGeom::uint;

void MeshLaplacian::constructTutte( const CMesh* tmesh )
{
	mOrder = tmesh->vertCount();
	std::vector<uint> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights(mOrder);

	for (int i = 0; i < mOrder; ++i) {
		const CVertex* vi = tmesh->getVertex(i);
		vector<int> vNeighbors;
		tmesh->vertRingNeighborVerts(i, 1, vNeighbors, false);
		int valence = vNeighbors.size();

		for (int j = 0; j < valence; ++j) {
			vII.push_back(i+1);
			vJJ.push_back(vNeighbors[j]+1);
			vSS.push_back(1.0);
		}
		vII.push_back(i+1);
		vJJ.push_back(i+1);
		vSS.push_back(-valence);
		vWeights[i] = valence;
	}

	mLS.convertFromCOO(mOrder, mOrder, vII, vJJ, vSS);
	mW.convertFromDiagonal(vWeights);

	mConstructed = true;
	m_laplacianType = Tutte;
}

void MeshLaplacian::constructUmbrella( const CMesh* tmesh )
{
	this->constructTutte(tmesh);

	mW.setToIdentity(mOrder);
	m_laplacianType = Umbrella;
}

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

	mConstructed = true;
	m_laplacianType = CotFormula;
}

void MeshLaplacian::constructSymCot( const CMesh* tmesh )
{
	constructCotFormula(tmesh);
	mW.setToIdentity(mOrder);
	m_laplacianType = SymCot;
}

void MeshLaplacian::constructFromMesh3( const CMesh* tmesh, int ringT, double hPara1, double hPara2 )
{
	// curvature difference based
	mOrder = tmesh->vertCount();

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;

	vector<std::tuple<int,int,double> > vSparseElements;

	ringT = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);
	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals_const();
	
	for (int vi = 0; vi < mOrder; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFaceIdx(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace(vFaces[fi]);
			double face_area = pfi->computeArea();
			for (int k = 0; k < 3; ++k)
			{
				int vki = pfi->getVertexIndex(k);
				if (vki == vi) continue;
				const CVertex* pvk = pfi->getVertex(k);

//				double w1 = std::exp(-std::pow(tmesh->getGeodesic(vi, vki), 2) / hPara1);
				double w1 = std::exp(-(pvi->getPosition()-pvk->getPosition()).length2() / hPara1);
//				double w2 = std::exp(-std::pow(pvi->getMeanCurvature() - pvk->getMeanCurvature(), 2) );// / hPara2);
				double w2 = std::exp(-std::pow(dotProduct3D(vVertNormals[vi], pvi->getPosition() - pvk->getPosition()), 2) / hPara2);

				double svalue = w1 * w2;
				svalue *= face_area;

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

	mConstructed = true;
}

void MeshLaplacian::constructFromMesh4(const CMesh* tmesh, int ringT, double hPara1, double hPara2)
{
	// similar to bilateral filtering
	mOrder = tmesh->vertCount();

	std::vector<int> vII, vJJ;
	std::vector<double> vSS;
	std::vector<double> vWeights;

	vector<std::tuple<int,int,double> > vSparseElements;

	ringT = 1;
	hPara1 = std::pow(tmesh->getAvgEdgeLength(), 2);
//	hPara2 = std::pow(tmesh->getAvgEdgeLength(), 2);
	hPara2 = tmesh->getAvgEdgeLength();
	const std::vector<Vector3D>& vVertNormals = tmesh->getVertNormals_const();

	for (int vi = 0; vi < mOrder; ++vi)
	{
		const CVertex* pvi = tmesh->getVertex(vi); 
		vector<int> vFaces = tmesh->getVertexAdjacentFaceIdx(vi, ringT);
		for (int fi = 0; fi < vFaces.size(); ++fi)
		{
			const CFace* pfi = tmesh->getFace(vFaces[fi]);
			double face_area = pfi->computeArea();
			for (int k = 0; k < 3; ++k)
			{
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
	mConstructed = true;
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
			double face_area = pfi->computeArea();
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
	mConstructed = true;
}
