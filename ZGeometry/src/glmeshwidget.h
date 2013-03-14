#pragma once
#include <GL/glew.h>
#include <QGLWidget>
#include <ZMesh.h>
#include <vector>
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"
#include "DiffusionShapeMatcher.h"

class GLMeshWidget : public QGLWidget
{
	Q_OBJECT

public:
	bool m_bShowLegend;
	bool m_bShowFeatures;
	bool m_bShowSignature;
	bool m_bShowRefPoint;
	bool m_bShowFeatureMatching;
	bool m_bDrawRegistration;
	bool m_bDrawMatching;
	bool m_bShowCorrespondenceLine;
	int	 m_nMeshLevel;

	enum {QZ_MOVE, QZ_PICK, QZ_DRAG} editMode;

	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();
	void fieldView(const Vector3D &center, const Vector3D &bbox);
	void addMesh(DifferentialMeshProcessor* pMP, RenderSettings* pRS);
	void setShapeMatcher(DiffusionShapeMatcher* p) { pDSM = p; }
	void setPointSize(double s) { m_dFeatureSphereRadius = s; }
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
	void setupObject(const CQrot& qrot, const Vector3D& trans) const;
	void drawLegend(QPainter* painter);
	void drawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* color);
	void drawMeshExt(const DifferentialMeshProcessor* pPM, const RenderSettings* renderSettings) const;
	void drawMatching(const DiffusionShapeMatcher* shapeMatcher, const RenderSettings* rs1, const RenderSettings* rs2) const;
signals:
	void vertexPicked1(int pRef);
	void vertexPicked2(int pRef);
private:
	void drawGL();
	void setupViewport(int width, int height);
	bool glPick(int x, int y, Vector3D& _p, int obj = 0);

	std::vector<DifferentialMeshProcessor*> vpMP;
	std::vector<RenderSettings*> vpRS;
	DiffusionShapeMatcher* pDSM;

	CArcball		g_arcball;
	GLfloat			g_EyeZ;
	GLdouble		g_myNear;
	GLdouble		g_myFar;
	GLdouble		g_myAngle;
	int				g_startx, g_starty;

	double			m_dFeatureSphereRadius;
	double			m_dMeshPointSize;
};

