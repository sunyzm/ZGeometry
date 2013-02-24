#pragma once
#include <gl/GL.h>
#include <Geometry.h>

const GLfloat preset_colors[][4] = {{0.53, 0.70, 0.93, 1.0}, {0.99, 0.73, 0.62, 1.0}};

class RenderSettings
{
public:
	RenderSettings() : mesh_color(preset_colors[0]), displayType(Mesh), 
		               showFeatures(false), showRefPoint(false), 
					   showColorSignature(false), selected(false), 
					   glPolygonMode(GL_FILL), display_shift(0, 0, 0) {}

	enum {PointCloud, Wireframe, Mesh, None} displayType;
	unsigned int glPolygonMode;
	bool showFeatures;
	bool showRefPoint;
	bool showColorSignature;
	bool selected;
	const GLfloat* mesh_color;
	Vector3D display_shift;
};