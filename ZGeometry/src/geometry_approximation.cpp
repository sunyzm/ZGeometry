#include "geometry_approximation.h"
#include <ppl.h>
#include <metis.h>
#include <ZGeom/util.h>
#include "differential_geometry.h"
#include "global.h"
#define USE_SPAMS

using std::vector;
using ZGeom::VecNd;
using ZGeom::logic_assert;
using ZGeom::runtime_assert;
using ZGeom::SparseApproxMethod;
using ZGeom::Dictionary;
using ZGeom::SparseCoding;

using ZGeom::MetisMeshPartition;

void calSGWDict(const ZGeom::EigenSystem& mhb, int waveletScaleNum, ZGeom::Dictionary& dict)
{
	ZGeom::DenseMatrixd matSGW;
    CStopWatch timer;
    timer.startTimer();
	computeSGWMat(mhb, waveletScaleNum, matSGW);
    timer.stopTimer("Wavelet time only: ");

	// normalize atoms to have norm 1
	const int vertCount = matSGW.colCount();
	const int totalAtomCount = matSGW.rowCount();
	dict.resize(totalAtomCount, vertCount);
    Concurrency::parallel_for(0, totalAtomCount, [&](int i) {
    	dict[i] = matSGW.getRowVec(i);
		dict[i].normalize(2.);
    });
}

void calHKDict(const ZGeom::EigenSystem& es, double timescale, ZGeom::Dictionary& dict)
{
    ZGeom::DenseMatrixd matHK = ZGeom::calHeatKernelMatrix(es, timescale);
	int vertCount = matHK.colCount(), atomCount = matHK.rowCount();

	dict.resize(atomCount, vertCount);
	for (int i = 0; i < atomCount; ++i) {
		dict[i] = matHK.getRowVec(i);
		dict[i].normalize(2.);
	}
}

int computePursuitCompressionBasisNum(int n, int m, int k1, int k2, double ratio)
{
	assert(m > 0);
	int nBasis1, nBasis2;
	nBasis1 = (3*n*k1*ratio) / (3*k2 + std::ceil(std::log(m)/std::log(2)));
	nBasis2 = (3*n*k1*ratio - m) / (3*k2);
	int nBasis = std::max(nBasis1, nBasis2);
	return nBasis;
}

void computeDictionary(DictionaryType dictType, const ZGeom::EigenSystem& es, ZGeom::Dictionary& dict)
{
	int vertCount = es.eigVecSize();
	int eigVecCount = es.eigVecCount();
	dict.clear();

	switch (dictType)
	{
	case DT_UNIT:
		dict.resize(vertCount, vertCount);
		for (int i = 0; i < vertCount; ++i) {
			dict[i].resize(vertCount, 0);
			dict[i][i] = 1.;
		}
		break;

	case DT_Fourier:
		dict.resize(eigVecCount, vertCount);
		for (int i = 0; i < eigVecCount; ++i)
			dict[i] = es.getEigVec(i);
		break;

	case DT_FourierSpikes:
		dict.resize(eigVecCount + vertCount, vertCount);
		for (int i = 0; i < eigVecCount; ++i)
			dict[i] = es.getEigVec(i);
		for (int i = 0; i < vertCount; ++i) {
			dict[eigVecCount + i].resize(vertCount, 0);
			dict[eigVecCount + i][i] = 1.;
		}
		break;

	case DT_SGW1:
		calSGWDict(es, 1, dict);
		break;
    case DT_SGW2:
        calSGWDict(es, 2, dict);
        break;
	case DT_SGW3:
		calSGWDict(es, 3, dict);
		break;
	case DT_SGW4:
		calSGWDict(es, 4, dict);
		break;
	case DT_SGW5:
		calSGWDict(es, 5, dict);
		break;
	case DT_SGW3MHB:
		calSGWDict(es, 3, dict);
		break;
	case DT_SGW4MHB:
		calSGWDict(es, 4, dict);
		break;
	case DT_SGW5MHB:
		calSGWDict(es, 5, dict);
		break;
	}

	if (dictType == DT_SGW3MHB || dictType == DT_SGW4MHB || dictType == DT_SGW5MHB) {
		int atomCount = dict.atomCount();
		dict.expandTo(atomCount + eigVecCount);
		for (int i = 0; i < eigVecCount; ++i)
			dict[atomCount + i] = es.getEigVec(i);
	}
}

vector<SparseCoding> multiExtractFront(const vector<SparseCoding>& vCodings, int n) 
{
	vector<SparseCoding> vNewCodings;
	for (auto sc : vCodings) vNewCodings.push_back(sc.extractFront(n));
	return vNewCodings;
}

vector<SparseCoding> multiExtractAfter(const vector<SparseCoding>& vCodings, int n)
{
	vector<SparseCoding> vNewCodings;
	for (auto sc : vCodings) vNewCodings.push_back(sc.extractBack(n));
	return vNewCodings;
}

void multiDictSparseDecomposition(
	const MeshCoordinates& coordInput, 
	const vector<const Dictionary*>& vDicts, 
	const vector<int>& vNNZ, 
	vector<vector<SparseCoding> > &vFinalCodings)
{
	const int vertCount = coordInput.size();
	vector<VecNd> vInputCoords = coordInput.to3Vec();

	const Dictionary &dict1 = *vDicts[0], &dict2 = *vDicts[1];
	Dictionary dictAll;
	ZGeom::combineDictionary(dict1, dict2, dictAll);
	int dictSize1 = dict1.size(), dictSize2 = dict2.size();
	int nnz1 = vNNZ[0], nnz2 = vNNZ[1];
	vFinalCodings.resize(2);	
	vector<SparseCoding> &vCodingsDict1 = vFinalCodings[0], &vCodingsDict2 = vFinalCodings[1];
	vector<SparseCoding> vInitCodings;
	ZGeom::SparseApproximationOptions approxOpts;
	approxOpts.mApproxMethod = ZGeom::SA_SOMP;
	approxOpts.mCodingSize = nnz1 + nnz2;

	vector<VecNd> vCoordDict1, vCoordDict2, vCoordInitApprox;
	multiChannelSparseApproximate(vInputCoords, dictAll, vInitCodings, approxOpts);
	multiChannelSparseReconstruct(dictAll, vInitCodings, vCoordInitApprox);
	std::cout << "- S-OMP error: " << coordInput.difference(MeshCoordinates(vertCount, vCoordInitApprox)) << "\n";
	auto vInitCodingFront = multiExtractFront(vInitCodings, dictSize1);
	std::cout << "- initial #MCA1: " << vInitCodingFront[0].size() 
		      << ", #MCA2: " << multiExtractAfter(vInitCodings, dictSize1)[0].size() << "\n";
	multiChannelSparseReconstruct(dictAll, vInitCodingFront, vCoordDict1);
	
	int numIter = 10;
	for (int k = 0; k < numIter; ++k)
	{
		multiChannelSparseApproximate(vCoordDict1, dict1, vCodingsDict1, approxOpts.setCodingSize(nnz1));
		multiChannelSparseReconstruct(dict1, vCodingsDict1, vCoordDict1);
		MeshCoordinates coordTMP1(vertCount, vCoordDict1);
		vCoordDict2 = coordInput.substract(coordTMP1).to3Vec();
		multiChannelSparseApproximate(vCoordDict2, dict2, vCodingsDict2, approxOpts.setCodingSize(nnz2));
		multiChannelSparseReconstruct(dict2, vCodingsDict2, vCoordDict2);
		MeshCoordinates coordTMP2(vertCount, vCoordDict2);
		MeshCoordinates coordTMP3 = coordInput.substract(coordTMP2);
		std::cout << "- Iteration " << k << " error: " << coordTMP3.difference(coordTMP1)<< "\n";
		if (k < numIter - 1) vCoordDict1 = coordTMP3.to3Vec();		
	}
}

ZGeom::VecNd singleChannelSparseInpaint(const ZGeom::VecNd& vSignal, const std::vector<int>& vMask, const ZGeom::Dictionary& dict, ZGeom::SparseCoding& sc)
{
    int signalSize = vSignal.size();
    int atomCount = dict.size();
    std::vector<double> vecSignalMasked;
    for (int k = 0; k < signalSize; ++k) {
        if (vMask[k] == 1) vecSignalMasked.push_back(vSignal[k]);
    }
    ZGeom::VecNd vSignalMasked(vecSignalMasked);
    int maskedSignalSize = vSignalMasked.size();
    ZGeom::Dictionary dictMasked;
    dictMasked.resize(atomCount, maskedSignalSize);
    for (int i = 0; i < atomCount; ++i) {
        const ZGeom::VecNd& vAtom = dict[i];
        std::vector<double> vecAtomMasked;
        for (int k = 0; k < signalSize; ++k) {
            if (vMask[k] == 1) vecAtomMasked.push_back(vAtom[k]);
        }
        dictMasked[i] = VecNd(vecAtomMasked);
    }

    ZGeom::SparseApproximationOptions approxOpts;
    approxOpts.mCodingSize = 50;
    approxOpts.mApproxMethod = ZGeom::SA_OMP;
    vector<SparseCoding> vOMPCodings;
    vector<VecNd> vApproximatedCoords;

    ZGeom::multiChannelSparseApproximate(vector < VecNd > {vSignalMasked}, dictMasked, vOMPCodings, approxOpts);
    sc = vOMPCodings[0];
    return ZGeom::singleChannelSparseReconstruct(dict, sc);
}

void ShapeApprox::init( CMesh* mesh )
{
	mOriginalMesh = mesh;
}

void ShapeApprox::doSegmentation( int maxSize )
{
	logic_assert(mOriginalMesh != NULL, "Error: Mesh is empty!");
	const int originalVertCount = mOriginalMesh->vertCount();
	int nPart = (maxSize > 0) ? std::lround((double)originalVertCount/(double)maxSize) : 1;
    if (nPart < 1) nPart = 1;

	if (nPart == 1)	// no segmentation performed; just copy the original mesh
	{
		mPartIdx.resize(originalVertCount, 0);		
		mSubMeshApprox.resize(1);
        mSubMeshApprox[0] = new SubMeshApprox;
		mSubMeshApprox[0]->mSubMesh.cloneFrom(*mOriginalMesh, ".sub0");
		mSubMeshApprox[0]->mMappedIdx.resize(originalVertCount);
		for (int i = 0; i < originalVertCount; ++i) 
			mSubMeshApprox[0]->mMappedIdx[i] = i;
		mSubMeshApprox[0]->init();
	}
	else 
	{
		mPartIdx = MetisMeshPartition(mOriginalMesh, nPart);
        //mPartIdx = ZGeom::FiedlerMeshPartition(mOriginalMesh, nPart);

		mSubMeshApprox.resize(nPart);
		std::vector<CMesh*> vSubMeshes;
		std::vector<std::vector<int>*> vMappedIdx;
		for (int partIdx = 0; partIdx < nPart; ++partIdx) {
            mSubMeshApprox[partIdx] = new SubMeshApprox;
			vSubMeshes.push_back(&mSubMeshApprox[partIdx]->mSubMesh);
			vMappedIdx.push_back(&mSubMeshApprox[partIdx]->mMappedIdx);
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
			mSubMeshApprox[s]->init();
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
	double decomposeTime = time_call_sec([&](){
		for (auto& m : mSubMeshApprox) {
			m->prepareEigenSystem(lapType, eigenCount);
		}
	});

	std::cout << "** Eigendecomposition finished! Time: " << decomposeTime << "s\n";
}

void ShapeApprox::constructDictionaries( DictionaryType dictType )
{
	double dictConstructTime = time_call_sec([=](){
		for (auto& m : mSubMeshApprox)
			m->constructDict(dictType);
	});
	std::cout << "** Dictionary constructed!  Time: " << dictConstructTime << "s\n";
}

void ShapeApprox::findSparseRepresentationBySize(SparseApproxMethod codingMethod, int codingSize )
{
	ZGeom::logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not segmented!");

	for (auto& m : mSubMeshApprox) {
		m->doSparseCoding(codingMethod, codingSize);
	}
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByRatio( SparseApproxMethod codingMethod, double basisRatio, bool encodeIndices )
{
	using ZGeom::logic_assert;
	logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	logic_assert(basisRatio > 0 && basisRatio <= 1., "!! Error: illegal coding ratio!");

	CStopWatch timer;
	timer.startTimer();
	int totalCodingSize(0);

	for (auto& m : mSubMeshApprox) {
		int codingSize = int(m->subMeshSize() * basisRatio);

		if (encodeIndices) {
			int subDictSize = m->dictSize();
			int indexBits = std::ceil(std::log((double)subDictSize)/std::log(2.));
			if (codingSize * indexBits < subDictSize) {
				codingSize += (subDictSize - codingSize*indexBits) / (96 + indexBits);
			}
		}

		if (codingSize > m->dictSize()) codingSize = m->dictSize();
		totalCodingSize += codingSize;
		m->doSparseCoding(codingMethod, codingSize);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByBasisRatio( SparseApproxMethod coding_method, double basisRatio )
{
	using ZGeom::logic_assert;
	logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	logic_assert(basisRatio > 0 && basisRatio <= 1., "!! Error: illegal coding ratio!");
	int totalCodingSize(0);
	CStopWatch timer;
	timer.startTimer();
	
	for (auto& m : mSubMeshApprox) {
		int coding_size = int(m->subMeshSize() * basisRatio);
		if (coding_size > m->dictSize()) coding_size = m->dictSize();
		totalCodingSize += coding_size;
		m->doSparseCoding(coding_method, coding_size);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::findSparseRepresentationByCompressionRatio( SparseApproxMethod codingMethod, double compressionRatio )
{
	using ZGeom::logic_assert;
	using ZGeom::runtime_assert;

	logic_assert(!mSubMeshApprox.empty(), "!! Error: Mesh is not partitioned!");
	logic_assert(compressionRatio > 0 && compressionRatio <= 1., "!! Error: illegal coding ratio!");
	int totalCodingSize(0);
	CStopWatch timer;
	timer.startTimer();

	for (auto& m : mSubMeshApprox) {
		int subDictSize = m->dictSize();
		int subMeshSize = m->subMeshSize();
		int k1 = 32, k2 = 32;	//bit size of float
		int coding_size = computePursuitCompressionBasisNum(subMeshSize, subDictSize, k1, k2, compressionRatio);
		runtime_assert(coding_size >= 1, "Number of participating basis must be positive!"); 

		if (coding_size > m->dictSize()) coding_size = m->dictSize();
		totalCodingSize += coding_size;
		m->doSparseCoding(codingMethod, coding_size);
	}

	timer.stopTimer("-- Encoding time: ");
	std::cout << "-- Total coding size: " << totalCodingSize << '\n';
	std::cout << "** Sparse Coding finished!\n";
}

void ShapeApprox::doSparseReconstructionBySize( int reconstructSize, MeshCoordinates& approxCoord )
{
	for (auto& m : mSubMeshApprox) {
		m->sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
}

void ShapeApprox::doSparseReconstructionByRatio( double basisRatio, MeshCoordinates& approxCoord, bool exploitSparsity )
{
	using ZGeom::logic_assert;

	logic_assert(0 < basisRatio && basisRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int reconstructSize = m->subMeshSize() * basisRatio;

		if (exploitSparsity) {
			int subDictSize = m->dictSize();
			int indexBits = int(std::log((double)subDictSize)/std::log(2.)) + 1;
			if (reconstructSize * indexBits < subDictSize) {
				reconstructSize += (subDictSize - reconstructSize*indexBits) / (96 + indexBits);
			}
		}

		if (reconstructSize > m->codingSize()) reconstructSize = m->codingSize();
		totalReconstructSize += reconstructSize;
		m->sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';
}

void ShapeApprox::doSparseReconstructionByBasisRatio( double basisRatio, MeshCoordinates& approxCoord )
{
	logic_assert(0 < basisRatio && basisRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int reconstruct_basis_num = m->subMeshSize() * basisRatio;
		if (reconstruct_basis_num > m->codingSize()) reconstruct_basis_num = m->codingSize();
		totalReconstructSize += reconstruct_basis_num;
		m->sparseReconstruct(reconstruct_basis_num);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';
}

void ShapeApprox::doSparseReconstructionByCompressionRatio( double compressionRatio, MeshCoordinates& approxCoord )
{
	logic_assert(0 < compressionRatio && compressionRatio <= 1.);
	int totalReconstructSize(0);

	for (auto& m : mSubMeshApprox) {
		int subDictSize = m->dictSize();
		int subMeshSize = m->subMeshSize();
		int k1 = 32, k2 = 32;	//bit size of float
		int reconstructSize = computePursuitCompressionBasisNum(subMeshSize, subDictSize, k1, k2, compressionRatio);
		runtime_assert(reconstructSize >= 1, "Number of participating basis must be positive!"); 

		if (reconstructSize > m->codingSize()) reconstructSize = m->codingSize();
		totalReconstructSize += reconstructSize;
		m->sparseReconstruct(reconstructSize);
	}
	integrateSubmeshApproximation(approxCoord);
	std::cout << "-- Total reconstruct size: " << totalReconstructSize << '\n';;
}

void ShapeApprox::doSparseReconstructionStepping( int totalSteps, std::vector<MeshCoordinates>& contCoords )
{
	contCoords.resize(totalSteps);

	for (int step = 0; step < totalSteps; ++step) {
		for (auto& m : mSubMeshApprox) m->sparseReconstructStep(step);
		integrateSubmeshApproximation(contCoords[step]);		 
	}
}

void ShapeApprox::integrateSubmeshApproximation(MeshCoordinates& integrated_approx_coord)
{
	/* integrate approximation results of all sub-meshes */
	const int vert_count = mOriginalMesh->vertCount();
	integrated_approx_coord.resize(vert_count);

	for (auto& m : mSubMeshApprox) {
		const std::vector<int>& vMappedIdx = m->mappedIdx();
		int sub_mesh_size = m->subMeshSize();

        for (int si = 0; si < sub_mesh_size; ++si) {
            int oi = vMappedIdx[si];
            for (int c = 0; c < 3; ++c)
                integrated_approx_coord(oi, c) = m->mReconstructedCoord(si, c);
        }
        
		//for (int c = 0; c < 3; ++c) {
		//	double* sub_approx_coord = m->mReconstructedCoord.getCoordData(c);
		//	for (int i = 0; i < sub_mesh_size; ++i) {
        //       integrated_approx_coord(vMappedIdx[i], c) = sub_approx_coord[i];
		//	}
		//}		
	}
}

void ShapeApprox::clear()
{
    mOriginalMesh = nullptr;
    mPartIdx.clear();
    //TODO: fix the following line
    //for (auto& sma : mSubMeshApprox) delete sma;
    mSubMeshApprox.clear();
}

void SubMeshApprox::prepareEigenSystem( LaplacianType lap_type, int eigen_num )
{
    MeshHelper& mMeshProcessor = mMeshHelper;
    mEigenSystem = mMeshHelper.prepareEigenSystem(lap_type, eigen_num);
}

void SubMeshApprox::constructDict( DictionaryType dictType )
{
	computeDictionary(dictType, mEigenSystem, mDict);
}

void SubMeshApprox::doSparseCoding( SparseApproxMethod approxMethod, int selected_atom_count )
{
	const int vertCount = mSubMesh.vertCount();
	const int atomCount = mDict.atomCount();
	runtime_assert(atomCount >= selected_atom_count);

    MeshCoordinates vertCoords = mSubMesh.getVertCoordinates();

    std::vector<ZGeom::VecNd> vSignals{ vertCoords.getXCoord(), vertCoords.getYCoord(), vertCoords.getZCoord() };
	ZGeom::SparseCoding vApproxX, vApproxY, vApproxZ;
    std::vector<ZGeom::SparseCoding*> vApproxCoeff{ &vApproxX, &vApproxY, &vApproxZ };
	
	for (int c = 0; c < 3; ++c) 
        mCoding[c].resize(selected_atom_count);

	if (approxMethod == ZGeom::SA_Truncation)
	{
		double innerProd[3];
		for (int i = 0; i < selected_atom_count; ++i) {
			for (int c = 0; c < 3; ++c)
				innerProd[c] = mDict[i].dot(vSignals[c]);
			for (int c = 0; c < 3; ++c)
				mCoding[c][i] = ZGeom::SparseCodingItem(i, innerProd[c]);
		}
	}
	else if (approxMethod == ZGeom::SA_SMP || approxMethod == ZGeom::SA_SOMP)
	{
		if (approxMethod == ZGeom::SA_SMP) {
			ZGeom::SimultaneousMP(vSignals, mDict.getAtoms(), selected_atom_count, vApproxCoeff);
		} else {
			ZGeom::SimultaneousOMP(vSignals, mDict.getAtoms(), selected_atom_count, vApproxCoeff);
		}

		for (int c = 0; c < 3; ++c) {
			for (int i = 0; i < selected_atom_count; ++i) {
				const ZGeom::SparseCodingItem& item = (*vApproxCoeff[c])[i];
				mCoding[c][i] = ZGeom::SparseCodingItem(item.index(), item.coeff());
			}
		}		
	}	
}

void SubMeshApprox::sparseReconstruct( int reconstruct_atom_num )
{
	int sub_vert_count = mSubMesh.vertCount();
	if (reconstruct_atom_num > this->codingSize()) reconstruct_atom_num = this->codingSize();
	if (reconstruct_atom_num <= 0) reconstruct_atom_num = this->codingSize();

	mReconstructedCoord.resize(sub_vert_count);
	for (int i = 0; i < reconstruct_atom_num; ++i) {
		for (int c = 0; c < 3; ++c) {
			const ZGeom::SparseCodingItem& sc = mCoding[c][i];
            const VecNd& selected_atom = mDict.getAtom(sc.index());
            for (int vi = 0; vi < sub_vert_count; ++vi)
                mReconstructedCoord(vi, c) += selected_atom[vi] * sc.coeff();
            
            // TODO: investigate the following lines
            //VecNd component_to_add = mDict.getAtom(sc.index());
            //component_to_add *= sc.coeff();            
            //mReconstructedCoord.addWith(c, component_to_add);
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
		const ZGeom::SparseCodingItem& sc = mCoding[c][step];
		mReconstructedCoord.addWith(c, sc.coeff() * mDict[sc.index()]);
	}
}

void SubMeshApprox::computeSparseCoding( const std::vector<double>& vecSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vCoeff )
{
	ZGeom::VecNd vSignal(vecSignal);
	computeSparseCoding(vSignal, opts, vCoeff);
}

void SubMeshApprox::computeSparseCoding( const ZGeom::VecNd& vSignal, SparseCodingOptions& opts, ZGeom::SparseCoding& vCoeff )
{
	int vertCount = mSubMesh.vertCount();
	int atomCount = mDict.atomCount();
	int codingAtomCount = opts.mCodingAtomCount;
	runtime_assert(atomCount >= codingAtomCount && vSignal.size() == vertCount);

	if (opts.mApproxMethod == ZGeom::SA_SMP || opts.mApproxMethod == ZGeom::SA_SOMP)
	{
		if (opts.mApproxMethod == ZGeom::SA_SMP) {
			ZGeom::MatchingPursuit(vSignal, mDict.getAtoms(), codingAtomCount, vCoeff);
		} else {
#ifdef USE_SPAMS
			ZGeom::OMP_SPAMS(g_engineWrapper, vSignal, mDict.getAtoms(), codingAtomCount, vCoeff);
#else
			ZGeom::OMP(vSignal, mDict.getAtoms(), codingAtomCount, vCoeff);
#endif
		}
	}
	else if (opts.mApproxMethod == ZGeom::SA_LASSO)
	{
		double lambda1 = opts.lambda1;
		ZGeom::LASSO_SPAMS(g_engineWrapper, vSignal, mDict.getAtoms(), lambda1, vCoeff);
	}
}
