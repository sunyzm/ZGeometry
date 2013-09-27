#pragma once
#include <iostream>
#include <ZMesh/Mesh.h>
#include <ZGeom/MatlabEngineWrapper.h>
#include "DifferentialMeshProcessor.h"

class ShapeEditor
{
public:
	ShapeEditor() : mMesh(nullptr), mProcessor(nullptr) {}

	void init(DifferentialMeshProcessor* processor) 
	{ 
		mProcessor = processor; 
		mMesh = processor->getMesh();
		processor->getMesh()->getVertCoordinates(mOldCoord); 
		mEngine = processor->getMatlabEngineWrapper();
		std::cout << "Shape editor is initialized!" << std::endl;
	}

	void revert() { mMesh->setVertCoordinates(mOldCoord); }

	void manifoldHarmonicsReconstruct(int nEig);
	void deformDifferential();
	void deformSimple();

	void deformSpectralWavelet();
	
private:
	void prepareAnchors(int& anchorCount, std::vector<int>& anchorIndex, std::vector<Vector3D>& anchorPos) const;

	CMesh* mMesh;
	MeshCoordinates mOldCoord;
	DifferentialMeshProcessor* mProcessor;
	ZGeom::MatlabEngineWrapper* mEngine;
	
};