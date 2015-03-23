#include <GL/glew.h>		// must include glew.h first
#include "glmeshwidget.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <QFile>
#include <QPainter>
#include <ZGeom/util.h>
#include <ZGeom/ZGeom.h>
#include <ZGeom/Color.h>
#include "OutputHelper.h"
#include "global.h"

using std::vector;
using ZGeom::Vec3d;
using ZGeom::Colorf;

extern OutputHelper qout;
extern GeometryTask g_task;

const GLfloat color_handle[4] = {1.0, 0.5, 0.5, 1};
const GLfloat featureColors[][4] = {{1.0, 0.0, 0.0, 1.0},	//red
									{0.0, 1.0, 0.0, 1.0},	//green
									{0.0, 0.0, 1.0, 0.0},	//blue
									{1.0, 1.0, 0.0, 1.0},	//yellow
									{1.0, 0.0, 1.0, 1.0},	//magenta
									{0.0, 1.0, 1.0, 1.0}};	//cyan
const int gFeatureColorNum = 6;
Qt::MouseButton gButton;

void glColorRed() { glColor4f(1.0, 0.0, 0.0, 1.0);  }
void glColorGreen() { glColor4f(0.0, 1.0, 0.0, 1.0); }
void glColorBlue() { glColor4f(0.0, 0.0, 1.0, 1.0); }
void glColorf(ZGeom::Colorf c) { glColor4f(c[0], c[1], c[2], c[3]);  }

void glColorCoded(float v, float pf)
{
	int ic = v;
	float f = v - ic;
	switch(ic)
	{
	case 0: glColor4f(1,f,0,pf); break;     // red -> yellow
	case 1: glColor4f(1-f,1,0,pf); break;   // yellow -> green
	case 2: glColor4f(0,1,f,pf); break;		// green -> cyan
	case 3: glColor4f(0,1-f,1,pf); break;	// cyan -> blue
	case 4: glColor4f(f,0,1,pf); break;     // blue -> purple 
	case 5: glColor4f(1,0,1-f,pf); break;   // purple -> red
	}
}

GLMeshWidget::GLMeshWidget(QWidget *parent) : QOpenGLWidget(parent)
{
	reset();
}

GLMeshWidget::~GLMeshWidget()
{
}

void GLMeshWidget::reset()
{
	g_EyeZ = 10.0;
	g_myNear = 1.0;
	g_myFar = 100.0;
	g_myAngle = 40.0;

	mBaseFeatureRadius = 0.02;
	mFeatureSphereRadius = mBaseFeatureRadius;
	mMeshPointSize = 1;
	m_nMeshLevel = 0;

	m_bShowLegend = false;
	m_bShowFeatures = false;
	m_bShowSignature = false;
	m_bShowLines = false;
	m_bShowRefPoint = false;
	m_bDrawMatching = false;
	m_bShowCorrespondenceLine = true;
	m_bDrawRegistration = false;
	m_bShowWireframeOverlay = false;
	m_bShowBoundingBox = false;
    m_bShowHoles = true;
    m_bShowSurrounding = true;
    m_bShowHoleError = true;
	
    m_nShadeMode = 1;
//	setAutoFillBackground(false);
}

void GLMeshWidget::initializeGL()
{
	qout.output("********************", OUT_TERMINAL);
    /* initialize GLEW */
    GLenum err = glewInit();
    if (GLEW_OK != err) {
        /* Problem: glewInit failed, something is seriously wrong. */
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

	/* print out some info about the graphics drivers */
	qout.output("OpenGL version: " + std::string((char *)glGetString(GL_VERSION)), OUT_TERMINAL);
	qout.output("GLSL version: " + std::string((char*)glGetString(GL_SHADING_LANGUAGE_VERSION)), OUT_TERMINAL);
	qout.output("Vendor: " + std::string((char*)glGetString(GL_VENDOR)), OUT_TERMINAL);
	qout.output("Renderer: " + std::string((char*)glGetString(GL_RENDERER)), OUT_TERMINAL);

    /* make sure OpenGL version 4.0 API is available */
	if(!GLEW_VERSION_4_0)
		qout.output("OpenGL 4.0 API is not available.", OUT_TERMINAL);
	qout.output("********************", OUT_TERMINAL);

    glEnable(GL_MULTISAMPLE);
}

void GLMeshWidget::resizeGL( int width, int height )
{
	setupViewport(width, height);
}

void GLMeshWidget::paintGL()
{
    drawGL();
}

void GLMeshWidget::mousePressEvent(QMouseEvent *event)
{
	const int win_width = this->width(), win_height = this->height();
	const int x = event->x(), y = event->y();

	std::vector<MeshHelper*>& vpMP = mMeshHelpers;
	std::vector<RenderSettings*>& vpRS = mRenderSettings;
	int meshCount = vpMP.size();

	if (editMode == QZ_MOVE)
	{
		if (event->button() == Qt::LeftButton) {
			gButton = Qt::LeftButton;
			g_arcball = CArcball(width(), height(), x - win_width/2, win_height/2 - y);
		} else if (event->button() == Qt::MidButton) {
			gButton = Qt::MidButton;
			g_startx = x;
			g_starty = y;
		} else if (event->button() == Qt::RightButton) {
			gButton = Qt::RightButton;
			g_startx = event->x();
			g_starty = event->y();
		}
	}
	
	else if (editMode == QZ_PICK)
	{
		if (event->button() == Qt::LeftButton) {
			/// for picking up a vertex
			Vec3d p;
			int obj_index = -1;

			if (meshCount >= 2) {
				if (vpRS[0]->selected && !vpRS[1]->selected) obj_index = 0;
				else if (!vpRS[0]->selected && vpRS[1]->selected) obj_index = 1;
			} 	else if (vpRS[0]->selected) obj_index = 0;

			if (obj_index >= 0 && glPick(x, y, p, obj_index)) {
				double dmin = 1e10;
				int hIdx = -1;
				int vertCount = vpMP[obj_index]->getMesh()->vertCount();
				for (int vi = 0; vi < vertCount; ++vi) {
					double d = (p - vpMP[obj_index]->getMesh()->vert(vi)->pos()).length();
					if (d < dmin) {
						dmin = d;
						hIdx = vi;
					}
				}
				qout.output("Pick coordinates: " + Int2String(x) + "," + Int2String(y), OUT_CONSOLE);
				qout.output("Pick vertex: #" + Int2String(hIdx), OUT_CONSOLE);

				if (event->modifiers() & Qt::ControlModifier) {
					if (obj_index == 0) {
						emit vertexPicked1(hIdx);	// change reference point
					} else if (obj_index == 1) {
						emit vertexPicked2(hIdx);
					}
				} else {
					vpMP[obj_index]->addNewHandle(hIdx);	// change handle set
				}				
			}
			// otherwise, no point picked up
		}
	}
	else if (editMode == QZ_DRAG)
	{
		/// for picking up a vertex
		int obj_index = -1;
		
		if (meshCount >= 2) {
			if (vpRS[0]->selected && !vpRS[1]->selected) obj_index = 0;
			else if (!vpRS[0]->selected && vpRS[1]->selected) obj_index = 1;
		} else if (vpRS[0]->selected) obj_index = 0;

		if (obj_index >= 0 && !vpMP[obj_index]->getHandles().empty()) {
			Vec3d p;
			if (glPick(x, y, p, obj_index)) {
				// find closest handle
				MeshHelper* pMP = vpMP[obj_index];
				int imin(-1);
				double d, dmin(1e10);
				for (auto handle : pMP->getHandles()) {
					d = handle.second.distanceTo(p);
					if (d < dmin) { 
						dmin = d; 
						imin = handle.first; 
					}
				}
				pMP->setActiveHandle(imin);
			}
		}
	}

	update();
}

void GLMeshWidget::mouseMoveEvent(QMouseEvent *event)
{
	const int win_width = this->width(), win_height = this->height();
	const int x = event->x(), y = event->y();
	std::vector<MeshHelper*>& vpMP = mMeshHelpers;
	std::vector<RenderSettings*>& vpRS = mRenderSettings;
	int meshCount = vpMP.size();

	if (editMode == QZ_PICK) return;
	else if (editMode == QZ_MOVE)
	{	
		ZGeom::Vec3d trans;
		CQrot    rot;
		
		if (gButton == Qt::LeftButton) {
			rot = g_arcball.update( x - win_width/2, win_height - y - win_height/2);
			for (int i = 0; i < meshCount; ++i) {
				if (vpRS[i]->selected)
					vpRS[i]->obj_rot = rot * vpRS[i]->obj_rot;
			}
		} else if (gButton == Qt::MidButton) {
			float scale = 3.0 * vpMP[0]->getMesh()->getBoundingBox().x / win_height;
			trans = ZGeom::Vec3d(scale * (x - g_startx), scale * (g_starty - y), 0);
			g_startx = x;
			g_starty = y;
			for (int i = 0; i < meshCount; ++i) {
				if (vpRS[i]->selected)
					vpRS[i]->obj_trans = vpRS[i]->obj_trans + trans;
			}
		}
		else if (gButton == Qt::RightButton ) {
			float scale = 5.0 * vpMP[0]->getMesh()->getBoundingBox().y / win_height;
			trans =  ZGeom::Vec3d(0, 0, scale * (g_starty - y));
			g_startx = x;
			g_starty = y;
			for (int i = 0; i < meshCount; ++i) {
				if (vpRS[i]->selected)
					vpRS[i]->obj_trans = vpRS[i]->obj_trans + trans;
			}
		}
	}
	else if (editMode == QZ_DRAG) 
	{
		int obj_index = -1;
		if (meshCount >= 2) {
			if (vpRS[0]->selected && !vpRS[1]->selected)
				obj_index = 0;
			else if (!vpRS[0]->selected && vpRS[1]->selected)
				obj_index = 1;
		}
		else if (meshCount == 1 && vpRS[0]->selected) 
			obj_index = 0;
		else return;

		if ( obj_index >= 0 && vpMP[obj_index]->getActiveHandle() != -1 
			 && vpMP[obj_index]->getHandles().find(vpMP[obj_index]->getActiveHandle()) != vpMP[obj_index]->getHandles().end() )
		{
			MeshHelper* pMP = vpMP[obj_index];

			GLdouble  modelview[16], projection[16];
			GLint     viewport[4];

			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			glLoadIdentity();
			gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
			CQrot& rot = vpRS[obj_index]->obj_rot;
			ZGeom::Vec3d& trans = vpRS[obj_index]->obj_trans;
			setupObject(rot, trans);

			glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
			glGetDoublev(GL_PROJECTION_MATRIX, projection);
			glGetIntegerv(GL_VIEWPORT, viewport);

			const ZGeom::Vec3d& handlePos = pMP->getHandles()[pMP->getActiveHandle()];
			int y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
			GLdouble ox, oy, oz, wx, wy, wz;
			gluProject(handlePos.x, handlePos.y, handlePos.z, modelview, projection, viewport, &wx, &wy, &wz);
			gluUnProject(x, y_new, wz, modelview, projection, viewport, &ox, &oy, &oz);
			pMP->getHandles()[pMP->getActiveHandle()] = ZGeom::Vec3d(ox, oy, oz);

			glPopMatrix();	
		}	
	}

	update();
}

void GLMeshWidget::mouseReleaseEvent( QMouseEvent *event )
{
}

void GLMeshWidget::wheelEvent(QWheelEvent *event)
{
	std::vector<RenderSettings*>& vpRS = mRenderSettings;
	int meshCount = mMeshHelpers.size();

	if (event->modifiers() & Qt::ControlModifier) {
		int numSteps = event->delta();
		float scale = 5.0 * mMeshHelpers[0]->getMesh()->getBoundingBox().x / this->height();
		ZGeom::Vec3d trans =  ZGeom::Vec3d(0, 0, scale * numSteps);
		
		for (int i = 0; i < meshCount; ++i) {
			if (vpRS[i]->selected)
				vpRS[i]->obj_trans = vpRS[i]->obj_trans + trans;
		}
	}

	update();
}

void GLMeshWidget::setupObject(const CQrot& qrot, const ZGeom::Vec3d& trans) const
{
	glMatrixMode(GL_MODELVIEW);
	glTranslated(trans.x, trans.y, trans.z);
	double rot[16];
	qrot.convert( rot );
	glMultMatrixd(( GLdouble*)rot );
}

void GLMeshWidget::fieldView( const ZGeom::Vec3d &center, const ZGeom::Vec3d &bbox )
{
	g_EyeZ = 15.0 * (float)bbox.y;
	g_myFar = 100.0 * (float)bbox.y;
	g_myNear = 0.01 * g_myFar;

	float len = bbox.x;
	if (bbox.y > len) len = bbox.y;
	if (bbox.z > len) len = bbox.z;

	g_myAngle = 2.0 * atan2(len, g_EyeZ);
	g_myAngle = (g_myAngle * 180.0) / ZGeom::PI + 2.0;
}

void GLMeshWidget::setupViewport( int width, int height )
{
	GLdouble ar = GLdouble(width) / GLdouble(height);	//aspect ratio

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(g_myAngle, ar, g_myNear, g_myFar);

//	glFrustum(-ar, ar, -1.0, 1.0, 4.0, 15.0);
// 	GLdouble clipX = g_myNear * tan(g_myAngle/2.0/180.0 * PI), 
// 		     clipY = clipX / ar;
// 	glFrustum(-clipX, -clipY, clipX, clipY, g_myNear, g_myFar);

	glMatrixMode(GL_MODELVIEW);
}


void GLMeshWidget::drawGL()
{
	static GLfloat position[] = {.0, .0, 1, 0.0};
	//static GLfloat diffuse[] = {1, 1, 1, 1};
	//static GLfloat global_ambient[] = {.2, .2, .2, 1};
	//static GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};

	makeCurrent();
    glClearColor(1., 1., 1., 0.);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushAttrib(GL_ALL_ATTRIB_BITS);
    	
    if (m_nShadeMode == 0) glShadeModel(GL_FLAT);
    else glShadeModel(GL_SMOOTH);

	glEnable(GL_DEPTH_TEST);
	glFrontFace(GL_CCW);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);		
	//glEnable (GL_BLEND); 
	//glBlendFunc (GL_SRC_ALPHA, GLblender_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
	//glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glLightfv(GL_LIGHT0, GL_POSITION, position);
	//glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	//glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	//glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);

	if (g_task == TASK_EDITING || g_task == TASK_VIEWING) {
        for (unsigned obj = 0; obj < mMeshHelpers.size(); ++obj) {
            drawMeshExt(mMeshHelpers[obj], mRenderSettings[obj]);
        }
	}
	else if (g_task == TASK_REGISTRATION) {
		assert(mMeshHelpers.size() == 2);
        drawMeshExt(mMatcher->getMeshHelper(0, m_nMeshLevel), mRenderSettings[0]);
        drawMeshExt(mMatcher->getMeshHelper(1, m_nMeshLevel), mRenderSettings[1]);

		if (m_bDrawMatching || m_bDrawRegistration) {
			drawCorrespondences(mMatcher, mRenderSettings[0], mRenderSettings[1]);
		}
	}

    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();

	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);
	glPopAttrib();
}

void GLMeshWidget::drawMeshExt( const MeshHelper* pMP, const RenderSettings* pRS ) const
{
    using ZGeom::MeshRegion;
    if (pMP->getMesh() == nullptr) return;
	CMesh* tmesh = pMP->getMesh();
    const vector<Vec3d>& vVertNormals = ZGeom::getMeshVertNormals(*tmesh);
    const vector<Vec3d> vVertPos = tmesh->allVertPos();
    const vector<Colorf>& vVertColors = tmesh->getVertColors(m_bShowSignature ? pRS->mActiveColorSignatureName : CMesh::StrAttrColorSigDefault);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	setupObject(pRS->obj_rot, pRS->obj_trans + pRS->display_shift);
	glPointSize(mMeshPointSize);
    
    //////////////////////////////////////////////////////////////////////////
	/* Start of drawing mesh */
    glPolygonMode(GL_FRONT_AND_BACK, pRS->glPolygonMode);
    glEnable(GL_POLYGON_OFFSET_FILL);

    /* draw mesh with color signature */
    {
        glPolygonOffset(1.0, 1.0);
        glBegin(GL_TRIANGLES);
        for (CFace* face : tmesh->m_vFaces) {
            for (int j = 0; j < 3; j++) {
                int pi = face->vertIdx(j);
                const ZGeom::Vec3d& norm = vVertNormals[pi];
                const Vec3d& vt = vVertPos[pi];
                const Colorf& vc = vVertColors[pi];
                glNormal3f(norm.x, norm.y, norm.z);
                glColor4f(vc[0], vc[1], vc[2], 1.0);
                glVertex3f(vt.x, vt.y, vt.z);
            }
        }
        glEnd();
    }

    /* draw hole regions in yellow */
    glPolygonOffset(0.5, 1.0);
    if (m_bShowHoles)
    {
        if (tmesh->hasAttr(StrAttrManualHoles))
        {
            vector<MeshRegion>& vManualHoles = tmesh->getAttrValue<vector<MeshRegion>>(StrAttrManualHoles);            
            for (MeshRegion& mr : vManualHoles) {
                const vector<int> &holeFaceIdx = mr.face_inside;                
                
                auto& vc = ZGeom::ColorYellow;
                glBegin(GL_TRIANGLES);                
                for (int fIdx : holeFaceIdx) {
                    CFace* face = tmesh->m_vFaces[fIdx];
                    for (int j = 0; j < 3; j++) {
                        int pi = face->vertIdx(j);
                        const ZGeom::Vec3d& norm = vVertNormals[pi];                        
                        const Vec3d& vt = vVertPos[pi];
                        glNormal3f(norm.x, norm.y, norm.z);
                        glColor4f(vc[0], vc[1], vc[2], 1.0f);
                        glVertex3f(vt.x, vt.y, vt.z);
                    }
                }
                glEnd();
            }   
        }
        else {
            vector<MeshRegion>& vHoles = tmesh->getAttrValue<vector<MeshRegion>>(ZGeom::StrAttrMeshHoleRegions);
            for (MeshRegion& mr : vHoles) {
                const vector<int> &holeFaceIdx = mr.face_inside;

                const vector<Colorf>* inpaint_error_colors = nullptr;
                if (m_bShowHoleError && tmesh->hasAttr(StrAttrColorInpaintError))
                    inpaint_error_colors = &tmesh->getVertColors(StrAttrColorInpaintError);

                glBegin(GL_TRIANGLES);                
                for (int fIdx : holeFaceIdx) {
                    CFace* face = tmesh->m_vFaces[fIdx];
                    for (int j = 0; j < 3; j++) {
                        int pi = face->vertIdx(j);
                        const ZGeom::Vec3d& norm = vVertNormals[pi];
                        Colorf vc = (inpaint_error_colors ? inpaint_error_colors->at(pi) : Colorf(ZGeom::ColorYellow));
                        const Vec3d& vt = vVertPos[pi];
                        glNormal3f(norm.x, norm.y, norm.z);                        
                        glColor4f(vc[0], vc[1], vc[2], 1.0f);
                        glVertex3f(vt.x, vt.y, vt.z);
                    }
                }
                glEnd();
            }
        }        
    }

    /* draw mesh surrounding faces */
    if (m_bShowSurrounding && tmesh->hasAttr(StrAttrHoleSurroundingFaces))
    {
        glPolygonOffset(0.5, 1.0);
        auto& surrounding_faces = tmesh->getAttrValue<vector<int>>(StrAttrHoleSurroundingFaces);
        glBegin(GL_TRIANGLES);
        const float *faceColor = ZGeom::ColorPink;
        glColor3f(faceColor[0], faceColor[1], faceColor[2]);
        for (int fIdx : surrounding_faces) {
            CFace* face = tmesh->m_vFaces[fIdx];
            for (int j = 0; j < 3; j++) {
                int pi = face->vertIdx(j);
                const ZGeom::Vec3d& norm = vVertNormals[pi];
                const Vec3d& vt = vVertPos[pi];
                const Colorf& vc = vVertColors[pi];
                glNormal3f(norm.x, norm.y, norm.z);
                glVertex3f(vt.x, vt.y, vt.z);
            }
        }
        glEnd();
    }

    glDisable(GL_POLYGON_OFFSET_FILL);   
    /* End of drawing mesh */
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    /* Start of drawing illustrative lines */
	glDisable(GL_LIGHTING); // disable lighting for overlaying lines
    glDisable(GL_LIGHT0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	/* highlight boundary edges */
#if 0
    if (tmesh->hasAttr(ZGeom::StrAttrMeshHoleRegions)) 
    {
        const vector<vector<int>>& boundaryLoops = ZGeom::getMeshBoundaryLoopHalfEdges(*tmesh);
        for (int i = 0; i < boundaryLoops.size(); ++i)
        {
            glBegin(GL_LINES);
            glLineWidth(4.0);
            if (boundaryLoops[i].size() <= ZGeom::MAX_HOLE_SIZE)
                glColor4f(0.0, 0, 1.0, 1.0);      // show edges on inner holes in blue
            else
                glColor4f(0.0, 0.0, 0.0, 1.0);	    // show edges on outer boundary in black
            for (int edgeIdx : boundaryLoops[i]) {
                const CHalfEdge* hf = tmesh->getHalfEdge(edgeIdx);
                int p1 = hf->getVertIndex(0), p2 = hf->getVertIndex(1);
                const Vec3d &v1 = vVertPos[p1], &v2 = vVertPos[p2];
                glVertex3d(v1.x, v1.y, v1.z);
                glVertex3d(v2.x, v2.y, v2.z);
            }
            glEnd();
        }
    }
#endif

    /* draw wireframe overlay */
    if (m_bShowWireframeOverlay)
    {
        auto color_overlay = ZGeom::ColorRed;
        glColor3f(color_overlay[0], color_overlay[1], color_overlay[2]);
        glPolygonOffset(1.0, 1.0);
        glBegin(GL_TRIANGLES);
        for (CFace* face : tmesh->m_vFaces) {
            for (int j = 0; j < 3; j++) {
                int pi = face->vertIdx(j);
                const ZGeom::Vec3d& norm = vVertNormals[pi];
                const Vec3d& vt = vVertPos[pi];
                glNormal3f(norm.x, norm.y, norm.z);
                glVertex3f(vt.x, vt.y, vt.z);
            }
        }
        glEnd();
    }

	/* draw vectors on mesh */
    if (m_bShowLines) 
    {
        for (std::string activeLineName : pRS->mActiveLineNames) {
            if (!tmesh->hasAttr(activeLineName)) continue;
            const MeshLineList& meshLines = tmesh->getAttrValue<MeshLineList>(activeLineName);
            const double avgEdgeLen = tmesh->getAvgEdgeLength();
            glLineWidth(1.0 * meshLines.lineWidthScale);
            glBegin(GL_LINES);
            for (auto line : meshLines) {
                if (line.directional) {
                    Vec3d v1 = line.first, vn = line.second;
                    Vec3d v2 = v1 + vn * avgEdgeLen * 0.5;
                    Vec3d vc = (v1 + v2) / 2.0;
                    glColorf(line.color1);	// vector line shooting from color1
                    glVertex3d(v1.x, v1.y, v1.z);
                    glVertex3d(vc.x, vc.y, vc.z);
                    glColorf(line.color2);	// to color2
                    glVertex3d(vc.x, vc.y, vc.z);
                    glVertex3d(v2.x, v2.y, v2.z);
                }
                else {
                    Vec3d v1 = line.first, v2 = line.second;
                    glColorf(line.color1);
                    glVertex3d(v1.x, v1.y, v1.z);
                    glVertex3d(v2.x, v2.y, v2.z);
                }
            }
            glEnd();
        }    
    }

    /* draw bounding box */
    if (m_bShowBoundingBox)
    {
        const ZGeom::Vec3d& bbox = tmesh->getBoundingBox();
        Vec3d verts[8] = { Vec3d(-bbox.x, -bbox.y, -bbox.z),
            Vec3d(-bbox.x, -bbox.y, bbox.z),
            Vec3d(-bbox.x, bbox.y, -bbox.z),
            Vec3d(-bbox.x, bbox.y, bbox.z),
            Vec3d(bbox.x, -bbox.y, -bbox.z),
            Vec3d(bbox.x, -bbox.y, bbox.z),
            Vec3d(bbox.x, bbox.y, -bbox.z),
            Vec3d(bbox.x, bbox.y, bbox.z) };
        auto drawQuad = [&verts](int i0, int i1, int i2, int i3) {
            glBegin(GL_LINE_LOOP);
            glVertex3d(verts[i0].x, verts[i0].y, verts[i0].z);
            glVertex3d(verts[i1].x, verts[i1].y, verts[i1].z);
            glVertex3d(verts[i2].x, verts[i2].y, verts[i2].z);
            glVertex3d(verts[i3].x, verts[i3].y, verts[i3].z);
            glEnd();
        };

        glLineWidth(1.0);
        glColor4f(0, 1.0, 1.0, 1.0);
        drawQuad(0, 2, 3, 1);
        drawQuad(1, 3, 7, 5);
        drawQuad(5, 7, 6, 4);
        drawQuad(4, 6, 2, 0);
        drawQuad(3, 2, 6, 7);
        drawQuad(0, 1, 5, 4);
    }

    glPolygonMode(GL_FRONT_AND_BACK, pRS->glPolygonMode);
    /* End of drawing illustrative lines */
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    /* Start of drawing illustrative points */
    glEnable(GL_LIGHTING);  
    glEnable(GL_LIGHT0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	/* draw reference point */
	if ( m_bShowRefPoint ) 
    {
        Vec3d vt = vVertPos[pMP->getRefPointIndex()];
		glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
		GLUquadric* quadric = gluNewQuadric();
		gluQuadricDrawStyle(quadric, GLU_FILL);
		glPushMatrix();
		glTranslated(vt.x, vt.y, vt.z);
		gluSphere(quadric, mFeatureSphereRadius, 16, 8);
		glPopMatrix();
	}

	/* draw feature points */
    if (m_bShowFeatures) 
    {
        for (std::string activePointFeature : pRS->mActivePointFeatures) {
            if (!tmesh->hasAttr(activePointFeature)) continue;
            const MeshFeatureList& feature_list = tmesh->getMeshFeatures(activePointFeature);
            /* draw as gluSphere */
            const float *feature_color1 = ZGeom::ColorGreen;
            const float *feature_color2 = ZGeom::ColorMagenta;
            GLUquadric* quadric = gluNewQuadric();
            for (MeshFeature* feature : feature_list.getFeatureVector()) {
                const Vec3d& vt = vVertPos[feature->m_index];
                glColor4f(feature->m_color[0], feature->m_color[1], feature->m_color[2], 1.0);                
                gluQuadricDrawStyle(quadric, GLU_FILL);
                glPushMatrix();
                glTranslated(vt.x, vt.y, vt.z);
                gluSphere(quadric, mFeatureSphereRadius, 16, 8);
                glPopMatrix();
            }
            gluDeleteQuadric(quadric);        
        }    
    } 

	/* draw handle points */
	if (!pMP->getHandles().empty())
	{
		for (auto handle :pMP->getHandles()) {
			ZGeom::Vec3d vt = handle.second;
			glColor4f(color_handle[0], color_handle[1], color_handle[2], color_handle[3]);
			GLUquadric* quadric = gluNewQuadric();
			gluQuadricDrawStyle(quadric, GLU_FILL);
			glPushMatrix();
			glTranslated(vt.x, vt.y, vt.z);
			gluSphere(quadric, mFeatureSphereRadius, 16, 8);
			glPopMatrix();
		}
	}

    glPolygonMode(GL_FRONT_AND_BACK, pRS->glPolygonMode);
    /* End of drawing illustrative points */
    //////////////////////////////////////////////////////////////////////////

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void GLMeshWidget::drawLegend(QPainter* painter)
{
	painter->setRenderHint(QPainter::Antialiasing);
	int xBegin = width() / 2 - 128;

	if (mLegendColors.empty()) return;
	for (int i = 0; i <= 255; i++)
	{
		//float gray = float(i) / float(255);
		//QColor col(255*FalseColorMap::red(gray), 255*FalseColorMap::green(gray), 255*FalseColorMap::blue(gray), 255);
		QColor col(255*mLegendColors[i].r(), 255*mLegendColors[i].g(), 255*mLegendColors[i].b(), 255);

		painter->setPen(QPen(col, 1, Qt::SolidLine));
		painter->drawLine(QPointF(xBegin+i, height()-50), QPointF(xBegin+i, height()-25));
	}
	painter->setPen(QPen(Qt::black, Qt::SolidLine));
//	painter->drawText(xBegin, height() - 70, 128, 12, Qt::AlignLeft, QString::number(mRenderSettings->at(0)->sigMin));
//	painter->drawText(xBegin + 128, height()-70, 128, 12, Qt::AlignRight, QString::number(mRenderSettings->at(0)->sigMax));
}

bool GLMeshWidget::glPick(int x, int y, ZGeom::Vec3d& _p, int obj /*= 0*/)
{
	std::vector<RenderSettings*>& vpRS = mRenderSettings;
	int meshCount = mMeshHelpers.size();
	if (obj >= meshCount || vpRS[obj]->selected == false) return false;

	GLdouble  modelview[16], projection[16];
	GLint     viewport[4];

    makeCurrent();
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
	const CQrot& rot = vpRS[obj]->obj_rot;
	const ZGeom::Vec3d& trans = vpRS[obj]->obj_trans;
	setupObject(rot, trans);

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

	// read depth buffer value at (x, y_new)
	float  z;
	int    y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
	glReadPixels(x, y_new, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );

	// reverse projection to get 3D point
	double pos[3];
	gluUnProject(x, y_new, z, modelview, projection, viewport, &pos[0], &pos[1], &pos[2]);

	glPopMatrix();

	if (z != 1.0f) {
		_p = ZGeom::Vec3d(pos[0], pos[1], pos[2]);
		return true;
	}
    else return false;
}

void GLMeshWidget::drawCorrespondences( const ShapeMatcher* shapeMatcher, const RenderSettings* rs1, const RenderSettings* rs2 ) const
{
	if (shapeMatcher == NULL || g_task != TASK_REGISTRATION) return;

	CMesh *tmesh1 = shapeMatcher->getMesh(0, 0), *tmesh2 = shapeMatcher->getMesh(1, 0);

	float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	glMateriali(GL_FRONT, GL_SHININESS, 96);
	double rot1[16], rot2[16];
	rs1->obj_rot.convert(rot1);
	rs2->obj_rot.convert(rot2);
	const ZGeom::Vec3d& trans1 = rs1->obj_trans;
	const ZGeom::Vec3d& trans2 = rs2->obj_trans;
		
	if (m_bDrawMatching)
	{
		const std::vector<MatchPair>& vmp = shapeMatcher->getMatchedFeaturesResults(shapeMatcher->getAlreadyMatchedLevel());
		const int size = (int)vmp.size();

		glLineWidth(2.0);
		GLUquadric* quadric = gluNewQuadric();
		gluQuadricDrawStyle(quadric, GLU_FILL);
		for (int i = 0; i < size; i++)
		{
			float cc = (i*1.0f)/(size-1.0f);
			glColorCoded(cc*4.0f, 0.8);
			//greenCoded(cc,0.9);
			int loc1 = vmp[i].m_idx1;
			int loc2 = vmp[i].m_idx2;

			const ZGeom::Vec3d& pos1 = tmesh1->vertPos(loc1);
			const ZGeom::Vec3d& pos2 = tmesh2->vertPos(loc2);

			double x1 = rot1[0]*pos1.x + rot1[4]*pos1.y + rot1[8]*pos1.z + trans1.x;
			double y1 = rot1[1]*pos1.x + rot1[5]*pos1.y + rot1[9]*pos1.z + trans1.y;
			double z1 = rot1[2]*pos1.x + rot1[6]*pos1.y + rot1[10]*pos1.z + trans1.z;

			double x2 = rot2[0]*pos2.x + rot2[4]*pos2.y + rot2[8]*pos2.z + trans2.x;
			double y2 = rot2[1]*pos2.x + rot2[5]*pos2.y + rot2[9]*pos2.z + trans2.y;
			double z2 = rot2[2]*pos2.x + rot2[6]*pos2.y + rot2[10]*pos2.z + trans2.z;

			glPushMatrix();
			glTranslated(x1, y1, z1);
			gluSphere(quadric, mFeatureSphereRadius, 16, 8);
			glPopMatrix();

			glPushMatrix();
			glTranslated(x2, y2, z2);
			gluSphere(quadric, mFeatureSphereRadius, 16, 8);
			glPopMatrix();

			if (m_bShowCorrespondenceLine) {
				if (vmp[i].m_note == -1)
// 					glColor4f(1.0, 0, 0, 1.0);
// 				else
					glColor4f(0.0, 0.0, 0.0, 1.0);
				glDisable(GL_LIGHTING);
				glBegin(GL_LINES);
				glVertex3d(x1, y1, z1);
				glVertex3d(x2, y2, z2);
				glEnd();
				glEnable(GL_LIGHTING);
			}
		}
		gluDeleteQuadric(quadric);
	}

	if (m_bDrawRegistration)
	{
		const std::vector<MatchPair>& vdr = shapeMatcher->getRegistrationResults(shapeMatcher->getAlreadyRegisteredLevel());
		const int size = vdr.size();
		std::vector<int> colorClass1, colorClass2;
		colorClass1.resize(tmesh1->vertCount(), -1);
		colorClass2.resize(tmesh2->vertCount(), -1);
		
		GLUquadric* quadric = gluNewQuadric();
		gluQuadricDrawStyle(quadric, GLU_FILL);
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		setupObject(rs1->obj_rot, rs1->obj_trans);	
		for (int i = 0; i < size; i++)
		{
			const int loc1 = vdr[i].m_idx1;
			const ZGeom::Vec3d& pos1 = tmesh1->vertPos(loc1);
			int color = i;
			float cc = (color*1.0f) / ((float)size-1.0f);

			glColorCoded(cc*4.0f, 0.9);
			glPushMatrix();
			glTranslated(pos1.x, pos1.y, pos1.z);
			gluSphere(quadric, tmesh1->getAvgEdgeLength()*0.2, 16, 8);
			glPopMatrix();
		}
		glPopMatrix();

		glPushMatrix();
		setupObject(rs2->obj_rot, rs2->obj_trans);	
		for (int i = 0; i < size; i++)
		{
			const int loc2 = vdr[i].m_idx2;
			const ZGeom::Vec3d& pos2 = tmesh2->vertPos(loc2);
			int color = i;
			float cc = (color*1.0f) / ((float)size-1.0f);

			glColorCoded(cc*4.0f, 0.9);
			glPushMatrix();
			glTranslated(pos2.x, pos2.y, pos2.z);
			gluSphere(quadric, tmesh2->getAvgEdgeLength()*0.2, 16, 8);
			glPopMatrix();
		}
		glPopMatrix();
		
		gluDeleteQuadric(quadric);
	}
}

void GLMeshWidget::changeShadeMode()
{
    m_nShadeMode = 1 - m_nShadeMode;
    update();
}

QImage GLMeshWidget::getScreenShot()
{
    return grabFramebuffer();
}
