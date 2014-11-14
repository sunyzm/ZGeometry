#include "Laplacian.h"
#include "EigenCompute.h"
#include "MeshProcessing.h"
#include "arithmetic.h"
#include "zassert.h"

namespace ZGeom {
using namespace std;

void Laplacian::decompose( int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, bool generalized /*= true*/ ) const
{
	runtime_assert(ep->isOpened(), "Matlab engine not opened!");
	EigenCompute eigenCompute(ep);
	if (generalized) {
		std::cout << "Do generalized eigendecomposition!\n";
		eigenCompute.solveGenSym(mLS, mW, nEig, eigSys);
	} else {
		std::cout << "Do standard eigendecomposition!\n";
		eigenCompute.solveStdSym(mLS, nEig, eigSys);
	}
}

void Laplacian::decomposeGeneralized(int nEig, MatlabEngineWrapper* ep, EigenSystem& eigSys, const SparseMatrix<double>& matB) const
{
	runtime_assert(ep->isOpened(), "Matlab engine not opened!");
	EigenCompute eigenCompute(ep);
	std::cout << "Do generalized eigendecomposition!\n";
	eigenCompute.solveGenSym(mLS, matB, nEig, eigSys);
}

void Laplacian::computeSubLaplacian( const std::vector<int>& vSelected, Laplacian& subLaplacian ) const
{
	const int subSize = (int)vSelected.size();
	subLaplacian.mOrder = subSize;
	mW.computeSubMatrix(vSelected, subLaplacian.mW);
	mLS.computeSubMatrix(vSelected, subLaplacian.mLS);

    //TODO: 
	//subLaplacian.mLS.makeLaplacian();
}

ZGeom::SparseMatrix<double> Laplacian::getSparseMatrix() const
{
    SparseMatrix<double> result = getLS();
    vector<double> vWeight = getW().getDiagonal();

    for (MatElem<double> &elem : result.allElements()) {
        elem.mVal /= vWeight[elem.row() - 1];
    }
    return result;
}

void Laplacian::constructTutte(const CMesh* tmesh)
{
    // L = D^-1 * (A - D)
    mOrder = tmesh->vertCount();
    ConstructMeshMatrix(*tmesh, ZGeom::MM_GRAPH_LAPLACE, mLS);
    mLS.scale(-1.0);
    ConstructMeshMatrix(*tmesh, ZGeom::MM_DEGREE, mW);
    mSymmetric = false;
}

void Laplacian::constructUmbrella(const CMesh* tmesh)
{
    // L = A - D
    mOrder = tmesh->vertCount();
    ConstructMeshMatrix(*tmesh, ZGeom::MM_GRAPH_LAPLACE, mLS);
    mLS.scale(-1);
    mW.setToIdentity(mOrder);	// set vertex weight matrix to identity; the attained Laplacian becomes symmetric
    mSymmetric = true;
}

void Laplacian::constructGeometricUmbrella(const CMesh *tmesh)
{
    cerr << "constructGeometricUmbrella not defined yet! " << endl;
    exit(-1);
    mSymmetric = true;
}

void Laplacian::constructNormalizedUmbrella(const CMesh* tmesh)
{
    /* L = D^(-1/2) * (D-A) * D^(-1/2) */
    mOrder = tmesh->vertCount();
    ConstructMeshMatrix(*tmesh, ZGeom::MM_NORMALIZED_GRAPH_LAPLACE, mLS);
    mLS.scale(-1);
    mW.setToIdentity(mOrder);
    mSymmetric = true;
}

/* Construct negative discrete Laplace operator */
void Laplacian::constructCotFormula(const CMesh* tmesh)
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

            double len0 = e0->length();
            double len1 = e1->length();
            double len2 = e2->length();
            amix += ZGeom::calMixedTriArea(len0, len1, len2);

            double cot_a11(0), cot_a12(0), cot_a21(0), cot_a22(0);
            ZGeom::triangleCot(len0, len1, len2, cot_a11, cot_a12);

            if (e0->twinHalfEdge() != NULL) {
                const CHalfEdge* e10 = e0->twinHalfEdge();
                const CHalfEdge* e11 = e10->nextHalfEdge();
                const CHalfEdge* e12 = e11->nextHalfEdge();
                len1 = e11->length();
                len2 = e12->length();

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
                len0 = e20->length();
                len1 = e21->length();
                len2 = e22->length();
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
    mSymmetric = false;
}

void Laplacian::constructSymCot(const CMesh* tmesh)
{
    const int vertCount = tmesh->vertCount();
    constructCotFormula(tmesh);
    mW.setToIdentity(vertCount);
    mSymmetric = true;
}

}   // end of namespace