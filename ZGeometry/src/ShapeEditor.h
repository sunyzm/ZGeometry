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

	void revert() { mMesh->setVertCoordinates(mOldCoord); };

	void manifoldHarmonicsReconstruct(int nEig);
	void differentialDeform();

private:
	CMesh* mMesh;
	MeshCoordinates mOldCoord;
	DifferentialMeshProcessor* mProcessor;
	ZGeom::MatlabEngineWrapper* mEngine;
	
};