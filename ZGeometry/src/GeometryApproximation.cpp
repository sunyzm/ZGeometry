#include "GeometryApproximation.h"
#include <metis.h>
#include "global.h"
using ZGeom::VecNd;

std::vector<int> MetisMeshPartition(const CMesh* mesh, int nPart)
{
	int vertCount = mesh->vertCount();
	std::vector<int> vPart(vertCount);
	int ncon = 1;
	std::vector<int> xadj, adjncy;
	mesh->getGraphCSR(xadj, adjncy);
	int objval;
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_CONTIG] = 1;
	options[METIS_OPTION_NUMBERING] = 0;

	METIS_PartGraphKway(&vertCount, &ncon, &xadj[0], &adjncy[0], NULL, NULL, NULL, &nPart, NULL, NULL, NULL, &objval, &vPart[0]);	
	
	return vPart;
}

void CalculateSGWDict(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::Dictionary& dict)
{
	ZGeom::DenseMatrixd matSGW;
	DifferentialMeshProcessor::computeSGWMat2(mhb, waveletScaleNum, matSGW);
	
	const int vertCount = matSGW.colCount();
	const int totalAtomCount = matSGW.rowCount();

	dict.resize(totalAtomCount, vertCount);
	for (int i = 0; i < totalAtomCount; ++i) {
		ZGeom::VecNd newBasis = matSGW.getRowVec(i);
		newBasis.normalize(ZGeom::RegularProductFunc);
		dict[i] = newBasis;
	}
}

int computePursuitCompressionBasisNum(int n, int m, int k1, int k2, double ratio)
{
	assert(m > 0);
	int nBasis1, nBasis2;
	nBasis1 = (3*n*k1*ratio) / (3*k2 + std::ceil(std::log(m)/std::log(2)));
	nBasis2 = (3*n*k1*ratio - m) / (3*k2);
	int nBasis = max(nBasis1, nBasis2);
	return nBasis;
}

void ShapeApprox::init( CMesh* mesh )
{
	mOriginalMesh = mesh;
}

void ShapeApprox::doSegmentation( int maxSize )
{
	ZUtil::logic_assert(mOriginalMesh != NULL, "Error: Mesh is empty!");
	const int originalVertCount = mOriginalMesh->vertCount();
	int nPart = (maxSize > 0) ? (originalVertCount / maxSize + 1) : 1;

	if (nPart == 1)	// no segmentation performed; just copy the original mesh
	{
		mPartIdx.resize(originalVertCount, 0);		
		mSubMeshApprox.resize(1);
		mSubMeshApprox[0].mSubMesh.cloneFrom(*mOriginalMesh, ".sub0");
		mSubMeshApprox[0].mMappedIdx.resize(originalVertCount);
		for (int i = 0; i < originalVertCount; ++i) 
			mSubMeshApprox[0].mMappedIdx[i] = i;
		mSubMeshApprox[0].init();
	}
	else 
	{
		mPartIdx = MetisMeshPartition(mOriginalMesh, nPart);		
		mSubMeshApprox.resize(nPart);
		std::vector<CMesh*> vSubMeshes;
		std::vector<std::vector<int>*> vMappedIdx;
		for (int partIdx = 0; partIdx < nPart; ++partIdx) {
			vSubMeshes.push_back(&mSubMeshApprox[partIdx].mSubMesh);
			vMappedIdx.push_back(&mSubMeshApprox[partIdx].mMappedIdx);
			vMappedIdx.back()->clear();
		}
		for (int vIdx = 0; vIdx < originalVertCount; ++vIdx) {
			vMappedIdx[mPartIdx[vIdx]]->push_back(vIdx);
		}

		mOriginalMesh->partitionToSubMeshes(vMappedIdx, vSubMeshes);
		
		for (int s = 0; s < nPart; ++s) {
			QString subMeshName;
			subMeshName.sprintf("%s.sub%d", mOriginalMesh->getMeshName().c_str(), s+1);
			vSubMeshes[s]->setMeshName(subMeshName.toStdString());
			mSubMeshApprox[s].init();
		}

		std::cout << "-- number of partitions: " << nPart << '\n';
		std::cout << "-- sub-mesh sizes: ";
		for (auto v: vMappedIdx) std::cout << v->size() << ' ';
		std::cout << '\n';
	}
	
	std::cout << "** Segmentation finished!" << std::endl;
}

void ShapeApprox::doEigenDecomposition( LaplacianType lapType, int eigenCount )
{
	CStopWatch timer;
	timer.startTimer();

	for (auto& m : mSubMeshApprox) {
		m.prepareEigenSystem(lapType, eigenCount);
	}

	timer.stopTimer("-- decomposition time: ");
	std::cout << "** Eigendecomposition finished!\n";
}

void ShapeApprox::constructDictionaries( DictionaryType dictType )
{
	for (auto& m : mSubMeshApprox)
		m.constructDict(dictType);
	std::cout << "** Dictionary constructed!\n";
}

void ShapeApprox::findSparseRepresentationBySize(SparseApproxMethod codingMethod, int codingSize )
{
	ZUtil::logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not segmented!");

	for (auto& m : mSubMeshApprox) {
		m.doSparseCoding(codingMethod, codingSize);
	}
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByRatio( SparseApproxMethod codingMethod, double basisRatio, bool encodeIndices )
{
	ZUtil::logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	ZUtil::logic_assert(basisRatio > 0 && basisRatio <= 1., "!! Error: illegal coding ratio!");

	CStopWatch timer;
	timer.startTimer();
	int totalCodingSize(0);

	for (auto& m : mSubMeshApprox) {
		int codingSize = int(m.subMeshSize() * basisRatio);

		if (encodeIndices) {
			int subDictSize = m.dictSize();
			int indexBits = std::ceil(std::log((double)subDictSize)/std::log(2.));
			if (codingSize * indexBits < subDictSize) {
				codingSize += (subDictSize - codingSize*indexBits) / (96 + indexBits);
			}
		}

		if (codingSize > m.dictSize()) codingSize = m.dictSize();
		totalCodingSize += codingSize;
		m.doSparseCoding(codingMethod, codingSize);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByBasisRatio( SparseApproxMethod codingMethod, double basisRatio )
{
	ZUtil::logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	ZUtil::logic_assert(basisRatio > 0 && basisRatio <= 1., "!! Error: illegal coding ratio!");
	int totalCodingSize(0);
	CStopWatch timer;
	timer.startTimer();
	
	for (auto& m : mSubMeshApprox) {
		int codingSize = int(m.subMeshSize() * basisRatio);
		if (codingSize > m.dictSize()) codingSize = m.dictSize();
		totalCodingSize += codingSize;
		m.doSparseCoding(codingMethod, codingSize);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByCompressionRatio( SparseApproxMethod codingMethod, double compressionRatio )
{
	ZUtil::logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	ZUtil::logic_assert(compressionRatio > 0 && compressionRatio <= 1., "!! Error: illegal coding ratio!");
	int totalCodingSize(0);
	CStopWatch timer;
	timer.startTimer();

	for (auto& m : mSubMeshApprox) {
		int subDictSize = m.dictSize();
		int subMeshSize = m.subMeshSize();
		int k1 = 32, k2 = 32;	//bit size of float
		int codingSize = computePursuitCompressionBasisNum(subMeshSize, subDictSize, k1, k2, compressionRatio);
		ZUtil::runtime_assert(codingSize >= 1, "Number of participating basis must be positive!"); 

		if (codingSize > m.dictSize()) codingSize = m.dictSize();
		totalCodingSize += codingSize;
		m.doSparseCoding(codingMethod, codingSize);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::doSparseReconstructionBySize( int reconstructSize, MeshCoordinates& approxCoord )
{
	for (auto& m : mSubMeshApprox) {
		m.sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
}

void ShapeApprox::doSparseReconstructionByRatio( double basisRatio, MeshCoordinates& approxCoord, bool exploitSparsity )
{
	ZUtil::logic_assert(0 < basisRatio && basisRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int reconstructSize = m.subMeshSize() * basisRatio;

		if (exploitSparsity) {
			int subDictSize = m.dictSize();
			int indexBits = int(std::log((double)subDictSize)/std::log(2.)) + 1;
			if (reconstructSize * indexBits < subDictSize) {
				reconstructSize += (subDictSize - reconstructSize*indexBits) / (96 + indexBits);
			}
		}

		if (reconstructSize > m.codingSize()) reconstructSize = m.codingSize();
		totalReconstructSize += reconstructSize;
		m.sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';
}

void ShapeApprox::doSparseReconstructionByBasisRatio( double basisRatio, MeshCoordinates& approxCoord )
{
	ZUtil::logic_assert(0 < basisRatio && basisRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int reconstructSize = m.subMeshSize() * basisRatio;

		if (reconstructSize > m.codingSize()) reconstructSize = m.codingSize();
		totalReconstructSize += reconstructSize;
		m.sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';
}

void ShapeApprox::doSparseReconstructionByCompressionRatio( double compressionRatio, MeshCoordinates& approxCoord )
{
	ZUtil::logic_assert(0 < compressionRatio && compressionRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int subDictSize = m.dictSize();
		int subMeshSize = m.subMeshSize();
		int k1 = 32, k2 = 32;	//bit size of float
		int reconstructSize = computePursuitCompressionBasisNum(subMeshSize, subDictSize, k1, k2, compressionRatio);
		ZUtil::runtime_assert(reconstructSize >= 1, "Number of participating basis must be positive!"); 

		if (reconstructSize > m.codingSize()) reconstructSize = m.codingSize();
		totalReconstructSize += reconstructSize;
		m.sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';;
}

void ShapeApprox::doSparseReconstructionStepping( int totalSteps, std::vector<MeshCoordinates>& contCoords )
{
	contCoords.resize(totalSteps);

	for (int step = 0; step < totalSteps; ++step) {
		for (auto& m : mSubMeshApprox) m.sparseReconstructStep(step);
		integrateSubmeshApproximation(contCoords[step]);		 
	}
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
	if (eigenCount == -1) eigenCount = mSubMesh.vertCount() - 1;

	int useCache = gSettings.LOAD_MHB_CACHE;
	if (useCache != 0 && mMeshProcessor.isMHBCacheValid(pathMHB, eigenCount)) {
		mMeshProcessor.loadMHB(pathMHB, laplacianType);
	}
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

	switch (dictType)
	{
	case DT_Fourier:
		mDict.resize(eigVecCount, vertCount);
		for (int i = 0; i < eigVecCount; ++i)
			mDict[i] = mEigenSystem.getEigVec(i);
		break;

	case DT_FourierSpikes:
		mDict.resize(eigVecCount + vertCount, vertCount);
		for (int i = 0; i < eigVecCount; ++i)
			mDict[i] = mEigenSystem.getEigVec(i);
		for (int i = 0; i < vertCount; ++i) {
			mDict[eigVecCount + i].resize(vertCount, 0);
			mDict[eigVecCount + i][i] = 1.;
		}
		break;

	case DT_SGW3:
		CalculateSGWDict(mEigenSystem, 3, mDict);
		break;
	case DT_SGW4:
		CalculateSGWDict(mEigenSystem, 4, mDict);
		break;
	case DT_SGW5:
		CalculateSGWDict(mEigenSystem, 5, mDict);
		break;
	case DT_SGW3MHB:
		CalculateSGWDict(mEigenSystem, 3, mDict);
		break;
	case DT_SGW4MHB:
		CalculateSGWDict(mEigenSystem, 4, mDict);
		break;
	case DT_SGW5MHB:
		CalculateSGWDict(mEigenSystem, 5, mDict);
		break;
	}

	if (dictType == DT_SGW3MHB || dictType == DT_SGW4MHB || dictType == DT_SGW5MHB) {
		int atomCount = mDict.atomCount();
		mDict.expandTo(atomCount + eigVecCount);
		for (int i = 0; i < eigVecCount; ++i)
			mDict[atomCount + i] = mEigenSystem.getEigVec(i);
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

// 		std::ofstream ofcoding("output/coding2.txt");
// 		for (int i = 0; i < selectedAtomCount; ++i) {
// 			for (int c = 0; c < 3; ++c) {
// 				auto& sc = mCoding[c][i];
// 				ofcoding << '(' << sc.mIdx << ',' << sc.mCoeff << ") "; 
// 			}
// 			ofcoding << '\n';
// 		}
	}	
}

void SubMeshApprox::sparseReconstruct( int reconstructAtomCount )
{
	int vertCount = mSubMesh.vertCount();
	if (reconstructAtomCount > this->codingSize()) reconstructAtomCount = this->codingSize();
	if (reconstructAtomCount <= 0) reconstructAtomCount = this->codingSize();

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