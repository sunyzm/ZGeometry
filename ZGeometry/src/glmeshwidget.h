#pragma once
#include <GL/glew.h>
#include <QGLWidget>
#include <ZMesh.h>
#include <vector>
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"


class GLMeshWidget : public QGLWidget
{
	Q_OBJECT

public:
	std::vector<DifferentialMeshProcessor*> vpMP;
	std::vector<RenderSettings> vSettings;

	bool m_bShowLegend;
	bool m_bShowFeatures;
	bool m_bShowSignature;
	bool m_bShowRefPoint;
	enum {QZ_MOVE, QZ_PICK, QZ_DRAG} editMode;

	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();
	void fieldView(const Vector3D &center, const Vector3D &bbox);
	void addMesh(DifferentialMeshProcessor* pmp);
protected:
	void initializeGL();
	void resizeGL(int width, int height);
//	void paintGL();
	void paintEvent(QPaintEvent *event);
//	void showEvent(QShowEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void setupObject(const CQrot& qrot, const Vector3D& trans);
	void drawLegend(QPainter* painter);
	void drawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* color);
	void drawMeshExt(int obj);
	void drawMeshExt(const DifferentialMeshProcessor* pPM, const Vector3D& trans, 
					 const CQrot& rot, const RenderSettings* renderSettings);

signals:
	void vertexPicked(int pRef);

private:
	void drawGL();
	void setupViewport(int width, int height);
	bool glPick(int x, int y, Vector3D& _p);

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

