#ifndef GLMESHWIDGET_H
#define GLMESHWIDGET_H

#include <QGLWidget>
#include <mesh/arcball.h>
#include "MeshProcessor.h"
#include <vector>

struct DisplaySettings
{
	DisplaySettings() : displayType(Mesh), showFeatures(false), showRefPoint(false) {}
	enum {PointCloud, Mesh, Signature} displayType;
	bool showFeatures;
	bool showRefPoint;
};

class GLMeshWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();
	void fieldView(const Vector3D &center, const Vector3D &bbox);

	std::vector<MeshProcessor*> vpMP;
	std::vector<DisplaySettings> vSettings;
	
	void addMesh(MeshProcessor* pmp);
protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void setupObject(const CQrot& qrot, const Vector3D& trans);
	void drawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* color);
private:
	void draw();
	
	int objSelect;	// -1 means all objects

	CArcball		g_arcball;
	GLfloat			g_EyeZ;
	CQrot			ObjRot1;
	Vector3D		ObjTrans1;
	CQrot			ObjRot2;
	Vector3D		ObjTrans2;
	GLdouble		g_myNear;
	GLdouble		g_myFar;
	GLdouble		g_myAngle;
	
	int				g_startx, g_starty;
};

#endif // GLMESHWIDGET_H
