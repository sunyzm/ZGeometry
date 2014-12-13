#pragma once
#include <vector>
#include <QGLWidget>
#include <QImage>
#include <ZGeom/Mesh.h>
#include <ZGeom/arcball.h>
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
	bool m_bShowLines;
	bool m_bShowRefPoint;
	bool m_bShowFeatureMatching;
	bool m_bDrawRegistration;
	bool m_bDrawMatching;
	bool m_bShowCorrespondenceLine;
	bool m_bShowWireframeOverlay;
	bool m_bShowBoundingBox;

    int  m_nShadeMode;
	int	 m_nMeshLevel;
	std::vector<ZGeom::Colorf> mLegendColors;

	enum { QZ_MOVE, QZ_PICK, QZ_DRAG } editMode;

	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();

	void setup(std::vector<DifferentialMeshProcessor*>* processors, std::vector<RenderSettings*>* rs, ShapeMatcher* matcher, ShapeEditor* editor) 
	{
		mProcessors = processors;
		mRenderSettings = rs;
		mMatcher = matcher;
		mEditor = editor;
	}

	void zoomPointSize(double s) { mFeatureSphereRadius = mBaseFeatureRadius * s; }
	void setBasePointSize(double r) { mFeatureSphereRadius *= r / mBaseFeatureRadius; mBaseFeatureRadius = r; }
	void fieldView(const Vector3D &center, const Vector3D &bbox);
	void reset();
    QImage getScreenShot();

protected:
	void initializeGL();
	void resizeGL(int width, int height);
    void paintGL();
    //void paintEvent(QPaintEvent *event);
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void setupObject(const CQrot& qrot, const Vector3D& trans) const;
	void drawMeshExt(const DifferentialMeshProcessor* pPM, const RenderSettings* renderSettings) const;
	void drawLegend(QPainter* painter);
	void drawCorrespondences(const ShapeMatcher* shapeMatcher, const RenderSettings* rs1, const RenderSettings* rs2) const;

signals:
	void vertexPicked1(int pRef);
	void vertexPicked2(int pRef);

public slots:
    void changeShadeMode();

private:
	void drawGL();
	void setupViewport(int width, int height);
    bool glPick(int x, int y, ZGeom::Vec3d& _p, int obj = 0);

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
	double			mBaseFeatureRadius;
	double			mFeatureSphereRadius;
	double			mMeshPointSize;
};

