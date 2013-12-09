#pragma once
#include <ZGeom/ZGeom.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/Mesh.h>
#include "DifferentialMeshProcessor.h"

class ShapeEditor
{
	friend class QZGeometryWindow;
public:
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}
	void init(DifferentialMeshProcessor* processor);
	void revertCoordinates();
	void changeCoordinates();
	void continuousReconstruct(int selectedApprox, int atomCount);
	const MeshCoordinates& getOldMeshCoord() const { return mCoords[0]; }

	void addNoise(double phi);

	void manifoldHarmonicsReconstruct(int nEig);
	void meanCurvatureFlow(double tMultiplier, int nRepeat = 1);

	void deformSimple();
	void deformLaplacian();
	void deformLaplacian_v2();
	void deformBiLaplacian();
	void deformMixedLaplacian(double ks, double kb);
	void deformThinShell2(double ks, double kb);
	void deformSpectralWavelet();	

private:
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void reconstructSpectralWavelet();
	void evalReconstruct(const MeshCoordinates& newCoord) const;
	
	void reconstructionTest1();
	void reconstructionTest2();
	void editTest2();

	void updateEditBasis(const std::vector<ZGeom::VecNd>& vAtoms, const std::vector<int>& vSelectedIdx);
	void computeApproximations(const std::vector<ZGeom::VecNd>& vAtoms, 
		                       ZGeom::FunctionApproximation* vApproxCoeff[3], 
							   int nReconstruct, 
							   std::vector<MeshCoordinates>& continuousCoords, 
							   MeshCoordinates& finalCoord);

	CMesh* mMesh;	
	DifferentialMeshProcessor* mProcessor;
	ZGeom::MatlabEngineWrapper* mEngine;
	std::vector<ZGeom::VecNd> mEditBasis;	


	ZGeom::FunctionApproximation mApproxCoeff[3];
	std::vector<MeshCoordinates> mContReconstructCoords[3];
	std::vector<MeshCoordinates> mCoords;
	int mCoordSelect;
};