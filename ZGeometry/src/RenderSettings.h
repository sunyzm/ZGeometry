#pragma once
#include <vector>
#include <set>
#include <iostream>
#include <ZGeom/Vec3.h>
#include <ZGeom/Quat.h>

class RenderSettings
{    
public:
    RenderSettings();

	enum { PointCloud, Wireframe, Mesh, None } displayType;

	unsigned int glPolygonMode;	
	bool selected;
	ZGeom::Vec3d obj_trans;
	CQrot obj_rot;
	ZGeom::Vec3d display_shift;
	std::string mActiveColorSignatureName;	
	std::set<std::string> mActivePointFeatures;
    std::set<std::string> mActiveLineNames;
};

