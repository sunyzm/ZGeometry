#pragma once
#include <ZMesh/Mesh.h>
#include "DifferentialMeshProcessor.h"

class ShapeEditor
{
public:
	void init(DifferentialMeshProcessor* processor) { mProcessor = processor; }
	void copy(MeshCoordinates& coords) const;
	
private:
	DifferentialMeshProcessor* mProcessor;
};