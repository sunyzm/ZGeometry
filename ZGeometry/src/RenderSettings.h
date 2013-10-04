#pragma once
#include <vector>
#include <ZMesh/Geometry.h>
#include <ZMesh/Quat.h>
#include <iostream>

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

	RenderSettings() : mesh_color(preset_mesh_colors[0]), displayType(Mesh), 
					   showFeatures(false), showRefPoint(false), 
					   showColorSignature(false), selected(false), 
					   glPolygonMode(Z_FILL), display_shift(0, 0, 0),
					   obj_rot(1,0,0,0), obj_trans(0,0,0) {}
	
	unsigned int glPolygonMode;
	bool showFeatures;
	bool showRefPoint;
	bool showColorSignature;
	bool selected;

	Vector3D obj_trans;
	CQrot obj_rot;
	Vector3D display_shift;
	const float* mesh_color;
	
	std::vector<double> vOriginalSignature;
	std::vector<double> vDisplaySignature;
	double sigMin, sigMax;

	void normalizeSignatureFrom(const std::vector<double>& vFrom);
	void logNormalizeSignatureFrom(const std::vector<double>& vFrom);
	void bandCurveSignatureFrom(const std::vector<double>& vFrom, double lowend, double highend);
	
	void importSignature(const std::vector<double>& vFrom);
	void normalizeSignature();
	void logNormalizeSignature();
	void bandCurveSignature(double lowend, double highend);
};

