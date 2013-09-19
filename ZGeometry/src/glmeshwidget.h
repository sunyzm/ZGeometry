#pragma once
#include <vector>
#include <QGLWidget>
#include <ZMesh/ZMesh.h>
#include "DifferentialMeshProcessor.h"
#include "RenderSettings.h"
#include "ShapeMatcher.h"
#include "ShapeEditor.h"

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

	void setup(std::vector<DifferentialMeshProcessor*>* processors, std::vector<RenderSettings*>* rs, ShapeMatcher* matcher, ShapeEditor* editor) {
		mProcessors = processors;
		mRenderSettings = rs;
		mMatcher = matcher;
		mEditor = editor;
	}

	void setPointSize(double s) { m_dFeatureSphereRadius = s; }
	void fieldView(const Vector3D &center, const Vector3D &bbox);

protected:
	void initializeGL();
	void resizeGL(int width, int height);
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
	void drawCorrespondences(const ShapeMatcher* shapeMatcher, const RenderSettings* rs1, const RenderSettings* rs2) const;

signals:
	void vertexPicked1(int pRef);
	void vertexPicked2(int pRef);

private:
	void drawGL();
	void setupViewport(int width, int height);
	bool glPick(int x, int y, Vector3D& _p, int obj = 0);

	std::vector<DifferentialMeshProcessor*>* mProcessors;
	std::vector<RenderSettings*>* mRenderSettings;
	ShapeMatcher* mMatcher;
	ShapeEditor* mEditor;

	CArcball		g_arcball;
	GLfloat			g_EyeZ;
	GLdouble		g_myNear;
	GLdouble		g_myFar;
	GLdouble		g_myAngle;
	int				g_startx, g_starty;
	double			m_dFeatureSphereRadius;
	double			m_dMeshPointSize;
};

