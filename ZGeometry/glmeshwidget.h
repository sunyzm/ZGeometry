#pragma once
#include <QGLWidget>
#include <ZMesh.h>
#include <vector>
#include "MeshProcessor.h"

struct DisplaySettings
{
	DisplaySettings() : displayType(Mesh), showFeatures(false), showRefPoint(false), showColorSignature(false), selected(false), glPolygonMode(GL_FILL) {}

	enum {PointCloud, Wireframe, Mesh, None} displayType;
	GLenum glPolygonMode;
	bool showFeatures;
	bool showRefPoint;
	bool showColorSignature;
	bool selected;
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
	bool m_bShowLegend;	
	void addMesh(MeshProcessor* pmp);
protected:
	void initializeGL();
	void resizeGL(int width, int height);
//	void paintGL();
	void paintEvent(QPaintEvent *event);
//	void showEvent(QShowEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void setupObject(const CQrot& qrot, const Vector3D& trans);
	void drawLegend(QPainter* painter);
	void drawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* color);
	void drawMeshExt(int obj);
	
private:
	void drawGL();
	void setupViewport(int width, int height);

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

