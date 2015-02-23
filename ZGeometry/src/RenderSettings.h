#pragma once
#include <vector>
#include <set>
#include <iostream>
#include <ZGeom/Vec3.h>
#include <ZGeom/Quat.h>

const int Z_POINT = 0x1B00;
const int Z_LINE  = 0x1B01;
const int Z_FILL  = 0x1B02;

const float preset_mesh_colors[][4] = {{0.53, 0.70, 0.93, 1.0}, 
									   {0.99, 0.73, 0.62, 1.0},
									   {0.63, 0.78, 0.63, 1.0}};

class RenderSettings
{    
public:
	enum { PointCloud, Wireframe, Mesh, None } displayType;

	RenderSettings() : 
		displayType(Mesh), selected(false), 
        glPolygonMode(Z_FILL), display_shift(0, 0, 0),
		obj_rot(1,0,0,0), obj_trans(0,0,0),
		mActiveColorSignatureName(StrAttrColorPreset)
	{
		
	}
	
	unsigned int glPolygonMode;	
	bool selected;
	ZGeom::Vec3d obj_trans;
	CQrot obj_rot;
	ZGeom::Vec3d display_shift;
	std::string mActiveColorSignatureName;	
	std::set<std::string> mActivePointFeatures;
    std::set<std::string> mActiveLineNames;
};

