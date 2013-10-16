#pragma once
#include <ZMesh/Mesh.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include <ZGeom/ZGeom.h>
#include "DifferentialMeshProcessor.h"

class ShapeEditor
{
	friend class QZGeometryWindow;
public:
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}
	void init(DifferentialMeshProcessor* processor);
	void revert() { mMesh->setVertCoordinates(mOldCoord); }

	void addNoise(double phi);

	void manifoldHarmonicsReconstruct(int nEig);
	void meanCurvatureFlow(double tMultiplier, int nRepeat = 1);

	void deformSimple();
	void deformLaplacian();
	void deformBiLaplacian();
	void deformMixedLaplacian(double ks, double kb);
	void deformThinShell2(double ks, double kb);
	void deformSpectralWavelet();	

private:
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;
	void reconstructSpectralWavelet();

	CMesh* mMesh;
	MeshCoordinates mOldCoord;
	DifferentialMeshProcessor* mProcessor;
	ZGeom::MatlabEngineWrapper* mEngine;	
};