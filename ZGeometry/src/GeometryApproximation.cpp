#include "GeometryApproximation.h"
using ZGeom::VecNd;


void ShapeApprox::init( CMesh* mesh )
{
	mOriginalMesh = mesh;
}

void ShapeApprox::doSegmentation( int maxSize )
{
	ZUtil::logic_assert(mOriginalMesh != NULL, "Error: Mesh is empty!");

	if (maxSize <= 0)	// no segmentation performed; just copy the original mesh
	{
		mSubMeshApprox.resize(1);
		mSubMeshApprox[0].mSubMesh.cloneFrom(*mOriginalMesh, ".sub0");
		int originalVertCount = mOriginalMesh->vertCount();
		mSubMeshApprox[0].mMappedIdx.resize(originalVertCount);
		for (int i = 0; i < originalVertCount; ++i) 
			mSubMeshApprox[0].mMappedIdx[i] = i;
		mSubMeshApprox[0].init();
	}

	std::cout << "Shape Approximation - Segmentation finished!" << std::endl;
}

void ShapeApprox::doEigenDecomposition( int eigenCount )
{
	for (auto& m : mSubMeshApprox) {
		m.prepareEigenSystem(Umbrella, eigenCount);
	}
	std::cout << "Shape Approximation - Preparation finished!\n";
}

void ShapeApprox::findSparseRepresentation( DictionaryType dictType, SparseApproxMethod codingMethod, int codingSize )
{
	ZUtil::logic_assert(!mSubMeshApprox.empty(), "Error: Mesh is not segmented!");
	int vertCount = mOriginalMesh->vertCount();
	int segmentationCount = mSubMeshApprox.size();

	for (auto& m : mSubMeshApprox) {
		m.constructDict(dictType);
		m.doSparseCoding(codingMethod, codingSize);
	}

	std::cout << "Shape Approximation - Sparse Coding finished!\n";
}

void ShapeApprox::sparseReconstruction( int reconstructSize )
{
	for (auto& m : mSubMeshApprox) {
		m.sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(mApproxCoord);
}

void ShapeApprox::sparseReconstructionStepping( int totalSteps, std::vector<MeshCoordinates>& contCoords )
{
	contCoords.resize(totalSteps);

	for (int step = 0; step < totalSteps; ++step) {
		for (auto& m : mSubMeshApprox) m.sparseReconstructStep(step);
		integrateSubmeshApproximation(contCoords[step]);		 
	}

	mApproxCoord = contCoords.back();
}

void ShapeApprox::integrateSubmeshApproximation(MeshCoordinates& integratedApproxCoord)
{
	/* integrate approximation results of all sub-meshes */
	const int vertCount = mOriginalMesh->vertCount();
	integratedApproxCoord.resize(vertCount);

	for (auto& m : mSubMeshApprox) {
		const std::vector<int>& vMappedIdx = m.mappedIdx();
		int subMeshSize = m.subMeshSize();
		for (int c = 0; c < 3; ++c) {
			VecNd& jointApproxCoord = integratedApproxCoord.getCoordFunc(c);
			const VecNd& subApproxCoord = m.mReconstructedCoord.getCoordFunc(c);
			for (int i = 0; i < subMeshSize; ++i) {
				jointApproxCoord[vMappedIdx[i]] = subApproxCoord[i];
			}
		}		
	}
}

void SubMeshApprox::prepareEigenSystem( LaplacianType laplacianType, int eigenCount )
{
	mMeshProcessor.constructLaplacian(laplacianType);
	std::string pathMHB = mMeshProcessor.generateMHBPath("cache/", laplacianType);
	if (mMeshProcessor.isMHBCacheValid(pathMHB, eigenCount))
		mMeshProcessor.loadMHB(pathMHB, laplacianType);
	else {
		mMeshProcessor.decomposeLaplacian(eigenCount, laplacianType);
		mMeshProcessor.saveMHB(pathMHB, laplacianType);
	}
	mEigenSystem = mMeshProcessor.getMHB(laplacianType);
}

void SubMeshApprox::constructDict( DictionaryType dictType )
{
	int vertCount = mSubMesh.vertCount();
	int eigVecCount = mEigenSystem.eigVecCount();	
	mDict.clear();

	if (dictType == DT_Fourier)
	{
		mDict.resize(eigVecCount, vertCount);
		for (int i = 0; i < eigVecCount; ++i)
			mDict[i] = mEigenSystem.getEigVec(i);
	} else if (dictType == DT_SGW1) {

	}
}

void SubMeshApprox::doSparseCoding( SparseApproxMethod approxMethod, int selectedAtomCount )
{
	const int vertCount = mSubMesh.vertCount();
	const int atomCount = mDict.atomCount();
	ZUtil::runtime_assert(atomCount >= selectedAtomCount);

	MeshCoordinates vertCoords;
	mSubMesh.getVertCoordinates(vertCoords);
	std::vector<ZGeom::VecNd> vSignals;
	vSignals.push_back(vertCoords.getXCoord()); 
	vSignals.push_back(vertCoords.getYCoord());
	vSignals.push_back(vertCoords.getZCoord());

	ZGeom::FunctionApproximation vApproxX, vApproxY, vApproxZ;
	std::vector<ZGeom::FunctionApproximation*> vApproxCoeff;
	vApproxCoeff.push_back(&vApproxX); 
	vApproxCoeff.push_back(&vApproxY); 
	vApproxCoeff.push_back(&vApproxZ);

	for (int c = 0; c < 3; ++c) mCoding[c].resize(selectedAtomCount);

	if (approxMethod == SA_Truncation)
	{
		double innerProd[3];
		for (int i = 0; i < selectedAtomCount; ++i) {
			for (int c = 0; c < 3; ++c)
				innerProd[c] = mDict[i].dot(vSignals[c]);
			for (int c = 0; c < 3; ++c)
				mCoding[c][i] = SparseCoeff(i, innerProd[c]);
		}
	}
	else if (approxMethod == SA_SMP || approxMethod == SA_SOMP)
	{
		if (approxMethod == SA_SMP) {
			ZGeom::SimultaneousMP(vSignals, mDict.getAtoms(), selectedAtomCount, vApproxCoeff);
		} else {
			ZGeom::SimultaneousOMP(vSignals, mDict.getAtoms(), selectedAtomCount, vApproxCoeff);
		}

		for (int c = 0; c < 3; ++c) {
			for (int i = 0; i < selectedAtomCount; ++i) {
				const ZGeom::ApproxItem& item = (*vApproxCoeff[c])[i];
				mCoding[c][i] = SparseCoeff(item.index(), item.coeff());
			}
		}				
	}	
}

void SubMeshApprox::sparseReconstruct( int reconstructAtomCount )
{
	int vertCount = mSubMesh.vertCount();
	int codingCoeffCount = mCoding[0].size();
	if (reconstructAtomCount > codingCoeffCount) reconstructAtomCount = codingCoeffCount;

	mReconstructedCoord.resize(vertCount);
	for (int i = 0; i < reconstructAtomCount; ++i) {
		for (int c = 0; c < 3; ++c) {
			const SparseCoeff& sc = mCoding[c][i];
			mReconstructedCoord.getCoordFunc(c) += sc.mCoeff * mDict[sc.mIdx];
		}
	}
}

void SubMeshApprox::sparseReconstructStep( int step )
{
	const int vertCount = mSubMesh.vertCount();
	const int codingCoeffCount = mCoding[0].size();

	if (step >= codingCoeffCount) return;
	if (step == 0) mReconstructedCoord.resize(vertCount);

	for (int c = 0; c < 3; ++c) {
		const SparseCoeff& sc = mCoding[c][step];
		mReconstructedCoord.getCoordFunc(c) += sc.mCoeff * mDict[sc.mIdx];
	}
}