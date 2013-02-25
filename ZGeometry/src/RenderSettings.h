#pragma once
#include <Geometry.h>
#include <Quat.h>
#include <vector>

#define Z_POINT 0x1B00
#define Z_LINE 0x1B01
#define Z_FILL 0x1B02

const float preset_colors[][4] = {{0.53, 0.70, 0.93, 1.0}, {0.99, 0.73, 0.62, 1.0}};

class RenderSettings
{
public:
	RenderSettings() : mesh_color(preset_colors[0]), displayType(Mesh), 
		               showFeatures(false), showRefPoint(false), 
					   showColorSignature(false), selected(false), 
					   glPolygonMode(Z_FILL), display_shift(0, 0, 0),
					   obj_rot(1,0,0,0), obj_trans(0,0,0) {}

	enum {PointCloud, Wireframe, Mesh, None} displayType;
	unsigned int glPolygonMode;
	bool showFeatures;
	bool showRefPoint;
	bool showColorSignature;
	bool selected;

	Vector3D obj_trans;
	CQrot obj_rot;
	Vector3D display_shift;
	
	const float* mesh_color;
	std::vector<double> vDisplaySignature;
	double sigMin, sigMax;

	void normalizeSignatureFrom(const std::vector<double>& vFrom);
	void logNormalizeSignatureFrom(const std::vector<double>& vFrom);
	void bandCurveSignatureFrom(const std::vector<double>& vFrom, double lowend, double highend);
};

