#pragma once
#include <vector>
#include <iostream>
#include <ZGeom/Geometry.h>
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
	enum {PointCloud, Wireframe, Mesh, None} displayType;

	RenderSettings();
	
	unsigned int glPolygonMode;
	bool showFeatures;
	bool showRefPoint;
	bool showColorSignature;
	bool selected;

	Vector3D obj_trans;
	CQrot obj_rot;
	Vector3D display_shift;
	const float* mesh_color;
	std::string mColorSignatureName;
};

