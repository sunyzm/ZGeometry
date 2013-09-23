#pragma once
#include <iostream>
#include <ZMesh/Mesh.h>
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
		std::cout << "Shape editor is initialized!" << std::endl;
	}

	void revert() { mMesh->setVertCoordinates(mOldCoord); };

	void addNewHandle( int hIdx )
	{
		auto iter = mHandles.find(hIdx);
		if (iter != mHandles.end()) mHandles.erase(iter);
		else mHandles[hIdx] = mMesh->getVertex(hIdx)->getPosition();	 
	}

	void manifoldHarmonicsReconstruct(int nEig);
	void differentialReconstruct(int nEig);

private:
	CMesh* mMesh;
	DifferentialMeshProcessor* mProcessor;
	MeshCoordinates mOldCoord;
	std::map<int, Vector3D> mHandles;
	int mActiveHandle;
};