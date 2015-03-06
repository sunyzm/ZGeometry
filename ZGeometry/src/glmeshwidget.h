#pragma once
#include <vector>
#include <QOpenGLWidget>
#include <QImage>
#include <ZGeom/Mesh.h>
#include <ZGeom/arcball.h>
#include "MeshHelper.h"
#include "RenderSettings.h"
#include "ShapeMatcher.h"
#include "ShapeEditor.h"

class GLMeshWidget : public QOpenGLWidget
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
    bool m_bShowHoles;

    int  m_nShadeMode;
	int	 m_nMeshLevel;
	std::vector<ZGeom::Colorf> mLegendColors;

	enum { QZ_MOVE, QZ_PICK, QZ_DRAG } editMode;

	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();

    void setup(int obj, MeshHelper* processor, RenderSettings* rs)
    {
        mMeshHelpers[obj] = processor;
        mRenderSettings[obj] = rs;
    }

    void setup(const std::vector<MeshHelper*>& processors, std::vector<RenderSettings>& vrs, ShapeMatcher* matcher) 
	{
        mMeshHelpers.clear();
        for (MeshHelper* mh : processors) mMeshHelpers.push_back(mh);
        mRenderSettings.clear();
        for (RenderSettings& rs : vrs) mRenderSettings.push_back(&rs);
		mMatcher = matcher;
        setBasePointSize(mMeshHelpers[0]->getMesh()->getAvgEdgeLength() / 4.);
	}

	void zoomPointSize(double s) { mFeatureSphereRadius = mBaseFeatureRadius * s; }
	void setBasePointSize(double r) { mFeatureSphereRadius *= r / mBaseFeatureRadius; mBaseFeatureRadius = r; }
	void fieldView(const ZGeom::Vec3d &center, const ZGeom::Vec3d &bbox);
	void reset();
    QImage getScreenShot();

protected:
	void initializeGL();
	void resizeGL(int width, int height);
    void paintGL();
	void mousePressEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
	void setupObject(const CQrot& qrot, const ZGeom::Vec3d& trans) const;
	void drawMeshExt(const MeshHelper* pPM, const RenderSettings* renderSettings) const;
	void drawLegend(QPainter* painter);
	void drawCorrespondences(const ShapeMatcher* shapeMatcher, const RenderSettings* rs1, const RenderSettings* rs2) const;

signals:
	void vertexPicked1(int pRef);
	void vertexPicked2(int pRef);

public slots:
    void changeShadeMode();

private:
    GLMeshWidget(const GLMeshWidget &);

private:
	void drawGL();
	void setupViewport(int width, int height);
    bool glPick(int x, int y, ZGeom::Vec3d& _p, int obj = 0);

	std::vector<MeshHelper*> mMeshHelpers;
	std::vector<RenderSettings*> mRenderSettings;
	ShapeMatcher* mMatcher;

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

