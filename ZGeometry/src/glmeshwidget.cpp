#include <cstdlib>
#include "glmeshwidget.h"
#include "OutputHelper.h"
#include <fstream>
#include <ZUtil.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <vector>
#include <QFile>

const GLfloat color1[4] = {0.53, 0.70, 0.93, 1.0};
const GLfloat color2[4] = {0.99, 0.73, 0.62, 1.0}; //{0.63,0.78,0.63,1.0};
const GLfloat color_handle[4] = {1.0, 0.5, 0.5, 1};
const GLfloat featureColors[][4] = {{1.0, 0.0, 0.0, 1.0}, 
                                    {0.0, 1.0, 0.0, 1.0}, 
                                    {0.0, 0.0, 1.0, 0.0},
                                    {1.0, 1.0, 0.0, 1.0}, 
                                    {1.0, 0.0, 1.0, 1.0},
                                    {0.0, 1.0, 1.0, 1.0}};
const int gFeatureColorNum = 6;
extern OutputHelper qout;
extern QString qformat;
Qt::MouseButton gButton;
FalseColorMap falseColorMap;

void glFalseColor(float v, float p)
{
	int floor = v * 255.0;
	glColor4f(falseColorMap.RedMap[floor], falseColorMap.GreenMap[floor], falseColorMap.BlueMap[floor], p);
}

GLMeshWidget::GLMeshWidget(QWidget *parent) : QGLWidget(parent)
{
//	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer | QGL::Rgba));
	g_EyeZ = 10.0;
	ObjRot1 = ObjRot2 = CQrot(1,0,0,0);
	ObjTrans1 = ObjTrans2 = Vector3D(0,0,0);
	g_myNear = 1.0;
	g_myFar = 100.0;
	g_myAngle = 40.0;

	m_bShowLegend = false;
	m_bShowFeatures = false;
	m_bShowSignature = false;
	m_bShowRefPoint = false;

	vSettings.resize(2, RenderSettings());

	setAutoFillBackground(false);
}

GLMeshWidget::~GLMeshWidget()
{

}

void GLMeshWidget::addMesh(DifferentialMeshProcessor* pmp)
{
	vpMP.push_back(pmp);
	if (vpMP.size() == 1) vSettings[0].mesh_color = preset_colors[0];
	else if (vpMP.size() == 2) vSettings[1].mesh_color =preset_colors[1];
	
	if (vpMP.size() == 2) vSettings[1].display_shift = Vector3D(vpMP[1]->getMesh()->getBoundingBox().x / 2., 0, 0);

}

void GLMeshWidget::mousePressEvent(QMouseEvent *event)
{
	const int win_width = this->width(), win_height = this->height();
	const int x = event->x(), y = event->y();

	if (editMode == QZ_MOVE)
	{
		if (event->button() == Qt::LeftButton)
		{
			gButton = Qt::LeftButton;
			g_arcball = CArcball(width(), height(), x - win_width/2, win_height/2 - y);
		}
		else if (event->button() == Qt::MidButton)
		{
			gButton = Qt::MidButton;
			g_startx = x;
			g_starty = y;
		}
		else if (event->button() == Qt::RightButton)
		{
			gButton = Qt::RightButton;
			g_startx = event->x();
			g_starty = event->y();
		}
	}
	
	else if (editMode == QZ_PICK)
	{
		if (event->button() == Qt::LeftButton)
		{
			/// for picking up a vertex
			Vector3D p;

			if (glPick(x, y, p))
			{
				double dmin = 1e10;
				int hIdx = -1;
				for (int vi = 0; vi < this->vpMP[0]->getMesh()->getVerticesNum(); ++vi)
				{
					double d = p.distantFrom(this->vpMP[0]->getMesh()->getVertex_const(vi)->getPosition());
					if (d < dmin) 
					{
						dmin = d;
						hIdx = vi;
					}
				}
				qout.output("Pick coordinates: " + Int2String(x) + "," + Int2String(y), OUT_CONSOLE);
				qout.output("Pick vertex: #" + Int2String(hIdx), OUT_CONSOLE);

				if (event->modifiers() & Qt::ControlModifier)
				{
					emit vertexPicked(hIdx);	// change reference point
				}
				else 
				{
					vpMP[0]->addNewHandle(hIdx);	// change handle set
				}				
			}
			// else, no point picked up
		}
	}
	else if (editMode == QZ_DRAG)
	{
		if (!vpMP[0]->mHandles.empty())
		{
			Vector3D p;
			if (glPick(x, y, p))
			{
				// find closest handle
				int imin(-1);
				double d, dmin(1e10);
				for (auto iter = vpMP[0]->mHandles.begin(); iter != vpMP[0]->mHandles.end(); ++iter)
				{
					d = iter->second.distantFrom(p);
					if (d < dmin) { dmin = d; imin = iter->first; }
				}
				vpMP[0]->active_handle = imin;
			}
		}
	}

	update();
}

void GLMeshWidget::mouseMoveEvent(QMouseEvent *event)
{
	const int win_width = this->width(), win_height = this->height();
	const int x = event->x(), y = event->y();
	
	if (editMode == QZ_PICK) return;
	else if (editMode == QZ_MOVE)
	{	
		Vector3D trans;
		CQrot    rot;
		
		if (gButton == Qt::LeftButton) 
		{
			rot = g_arcball.update( x - win_width/2, win_height - y - win_height/2);
			if (vSettings[0].selected) 
				ObjRot1 = rot * ObjRot1;
			if (vSettings[1].selected) 
				ObjRot2 = rot * ObjRot2;
		}
		else if (gButton == Qt::MidButton) 
		{
			float scale = 3.0 * vpMP[0]->getMesh()->getBoundingBox().x / win_height;
			trans = Vector3D(scale * (x - g_startx), scale * (g_starty - y), 0);
			g_startx = x;
			g_starty = y;
			if (vSettings[0].selected) 
				ObjTrans1 = ObjTrans1 + trans;
			if (vSettings[1].selected) 
				ObjTrans2 = ObjTrans2 + trans;
		}
		else if (gButton == Qt::RightButton ) 
		{
			float scale = 5.0 * vpMP[0]->getMesh()->getBoundingBox().y / win_height;
			trans =  Vector3D(0, 0, scale * (g_starty - y));
			g_startx = x;
			g_starty = y;
			if (vSettings[0].selected)
				ObjTrans1 = ObjTrans1 + trans;
			if (vSettings[1].selected) 
				ObjTrans2 = ObjTrans2 + trans;
		}
	}
	else if (editMode == QZ_DRAG && vpMP[0]->active_handle != -1 && vpMP[0]->mHandles.find(vpMP[0]->active_handle) != vpMP[0]->mHandles.end())
	{
		
		GLdouble  modelview[16], projection[16];
		GLint     viewport[4];

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
		CQrot& rot = ObjRot1;
		Vector3D& trans = ObjTrans1;
		setupObject(rot, trans);

		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, projection);
		glGetIntegerv(GL_VIEWPORT, viewport);

		const Vector3D& handlePos = vpMP[0]->mHandles[vpMP[0]->active_handle];
		int y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
		GLdouble ox, oy, oz, wx, wy, wz;
		gluProject(handlePos.x, handlePos.y, handlePos.z, modelview, projection, viewport, &wx, &wy, &wz);
		gluUnProject(x, y_new, wz, modelview, projection, viewport, &ox, &oy, &oz);
		vpMP[0]->mHandles[vpMP[0]->active_handle] = Vector3D(ox, oy, oz);

		glPopMatrix();		
	}

	update();
}

void GLMeshWidget::mouseReleaseEvent( QMouseEvent *event )
{
	//if (editMode == QZ_DRAG)
	//{
	//	active_handle = -1;
	//}
}

void GLMeshWidget::wheelEvent(QWheelEvent *event)
{
	if (event->modifiers() & Qt::ControlModifier)
	{
		int numSteps = event->delta();
		float scale = 3.0 * vpMP[0]->getMesh()->getBoundingBox().x / this->height();
		Vector3D trans =  Vector3D(0, 0, scale * numSteps);
		
		if (vSettings[0].selected) 
			ObjTrans1 = ObjTrans1 + trans;
		if (vSettings[1].selected)
			ObjTrans2 = ObjTrans2 + trans;		
	}

	update();
}

void GLMeshWidget::initializeGL()
{
	// initialize GLEW
	qout.output("********************");
	if(glewInit() != GLEW_OK)
		qout.output("glewInit failed", OUT_CONSOLE);
	else qout.output("glewInit succeeded", OUT_CONSOLE);
	// print out some info about the graphics drivers

	qout.output("OpenGL version: " + std::string((char *)glGetString(GL_VERSION)), OUT_CONSOLE);
	qout.output("GLSL version: " + std::string((char*)glGetString(GL_SHADING_LANGUAGE_VERSION)), OUT_CONSOLE);
	qout.output("Vendor: " + std::string((char*)glGetString(GL_VENDOR)), OUT_CONSOLE);
	qout.output("Renderer: " + std::string((char*)glGetString(GL_RENDERER)), OUT_CONSOLE);

//#define MORE_DEBUG
#ifdef MORE_DEBUG	
	std::string strExt = std::string((char*)glGetString(GL_EXTENSIONS));
	std::ofstream ofs("output/glExt.txt");
	ofs << strExt;
	ofs.close();
//	std::system("python python/sp2ln.py output/glExt.txt");
#endif
	// make sure OpenGL version 3.2 API is available
	if(!GLEW_VERSION_3_2)
		qout.output("OpenGL 3.2 API is not available.", OUT_CONSOLE);
	qout.output("********************");
}

void GLMeshWidget::fieldView( const Vector3D &center, const Vector3D &bbox )
{
	g_EyeZ = 15.0 * (float)bbox.y;
	g_myFar = 100.0 * (float)bbox.y;
	g_myNear = 0.01 * g_myFar;

	float len = bbox.x;
	if (bbox.y > len) len = bbox.y;
	if (bbox.z > len) len = bbox.z;

	g_myAngle = 2.0 * atan2(len, g_EyeZ);
	g_myAngle = (g_myAngle * 180.0) / PI + 2.0;
}

void GLMeshWidget::resizeGL( int width, int height )
{
	setupViewport(width, height);
}

void GLMeshWidget::drawGL()
{
	makeCurrent();
 	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glMatrixMode(GL_MODELVIEW);
 	glPushMatrix();

	GLfloat diffuse[] = {1, 1, 1, 1};
	GLfloat global_ambient[] = {.2, .2, .2, 1};
	GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat position[] = {.0, .0, 1, 0.0};

	glClearColor(1., 1., 1., 0.);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glFrontFace(GL_CCW);
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
//	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
	glEnable (GL_BLEND); 
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_POSITION, position);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	setupViewport(width(), height());
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);

	if (vpMP.size() > 0 && vpMP[0] != NULL)
		drawMeshExt(vpMP[0], ObjTrans1, ObjRot1, &vSettings[0]);
	if (vpMP.size() > 1 && vpMP[1] != NULL)
		drawMeshExt(vpMP[1], ObjTrans2, ObjRot2, &vSettings[1]);

 	glMatrixMode(GL_MODELVIEW);
 	glPopMatrix();
 	glPopAttrib();
}

void GLMeshWidget::setupObject(const CQrot& qrot, const Vector3D& trans)
{
	//prerequisite: glMatrixMode(GL_MODELVIEW)
	glTranslated(trans.x, trans.y, trans.z);
	double rot[16];
	qrot.convert( rot );
	glMultMatrixd(( GLdouble*)rot );
}

/*
void GLMeshWidget::drawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* color)
{
	if(!tmesh) return;

	float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
//  glMateriali(GL_FRONT, GL_SHININESS, 96);

	glPushMatrix();
	setupObject(rot, trans);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	{	// just display mesh in single color
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tmesh->getFaceNum(); i++)
		{
			if(NULL == tmesh->m_pFace[i].m_piEdge) continue;
			for (int j = 0; j < 3; j++)
			{
				glColor4f(color[0], color[1], color[2], color[3]); 
				int pi = tmesh->m_pFace[i].m_piVertex[j];
				Vector3D norm = tmesh->m_pVertex[pi].getNormal();
				glNormal3f(norm.x, norm.y, norm.z);
				Vector3D vt = tmesh->m_pVertex[pi].getPosition();
				//vt -= tmesh->m_Center;
				glVertex3f(vt.x, vt.y, vt.z);
			}
		}
		glEnd();
	}
	glDisable(GL_POLYGON_OFFSET_FILL);


	if (tmesh->hasBounary())   //highlight boundary edge 
	{
		glDisable(GL_LIGHTING);
		glBegin(GL_LINES);	
		for(int i = 0; i < tmesh->getHalfEdgeNum(); i++)
		{
			if(tmesh->m_pHalfEdge[i].m_iTwinEdge < 0) 
			{
				int p1 = tmesh->m_pHalfEdge[i].m_iVertex[0];
				int p2 = tmesh->m_pHalfEdge[i].m_iVertex[1];
				glLineWidth(2.0);
				glColor4f(0.0, 0.0, 0.0, 1.0);			//show boundary edge in black
				if(tmesh->m_pVertex[p1].m_bIsHole) 
				{
					glColor4f(0.0, 0.0, 1.0, 1.0);		//show edge on holes in blue
				}
				Vector3D v1 = tmesh->getVertex_const(p1)->getPosition();
				//v1 -= tmesh->m_Center;
				Vector3D v2 = tmesh->getVertex_const(p2)->getPosition();
				//v2 -= tmesh->m_Center;
				glVertex3d(v1.x, v1.y, v1.z);
				glVertex3d(v2.x, v2.y, v2.z);
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
	
	glPopMatrix();
}
*/

#ifdef USE_OLD_DRAW_MESH
void GLMeshWidget::drawMeshExt( int obj )
{
	if (obj >= vpMP.size() || obj < 0) return;	
	if(!vpMP[obj]->getMesh()) return;

	const CMesh* tmesh = vpMP[obj]->getMesh();
	CQrot rot = (obj == 0) ? ObjRot1 : ObjRot2;
	Vector3D trans = (obj == 0) ? ObjTrans1 : ObjTrans2;
	const GLfloat *color = preset_colors[obj];	

	float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
//  glMateriali(GL_FRONT, GL_SHININESS, 96);

	Vector3D shift = Vector3D(0, 0, 0);
	if (obj == 1) shift = Vector3D(tmesh->getBoundingBox().x/2, 0, 0);

	bool showSignature = m_bShowSignature && !vpMP[obj]->vDisplaySignature.empty();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	setupObject(rot, trans);	

	/// draw the mesh
 	glPolygonMode(GL_FRONT_AND_BACK, vSettings[obj].glPolygonMode);
	
	glPointSize(2.0);
	glColor4f(color[0], color[1], color[2], color[3]); 
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);
	
	glBegin(GL_TRIANGLES);
	for (int i = 0; i < tmesh->getFaceNum(); i++)
	{
		if(!tmesh->m_pFace[i].hasHalfEdge()) continue;
		for (int j = 0; j < 3; j++)
		{
			int pi = tmesh->m_pFace[i].getVertexIndex(j);					 				
			Vector3D norm = tmesh->getVertex_const(pi)->getNormal();
			glNormal3f(norm.x, norm.y, norm.z);	 				
			Vector3D vt = tmesh->getVertex_const(pi)->getPosition();
			vt += shift;	//add some offset to separate object 1 and 2
			if (showSignature) 
				glFalseColor(vpMP[obj]->vDisplaySignature[pi], 1.0);
			glVertex3f(vt.x, vt.y, vt.z);
		}
	}
	glEnd();
	
	glDisable(GL_POLYGON_OFFSET_FILL);
	
#ifdef MORE_DRAWING
	/// draw boundary edge
//  	glDisable(GL_LIGHTING);
//   	if (tmesh->hasBounary())   //highlight boundary edge 
//   	{
//   		glBegin(GL_LINES);	
//   		for(int i = 0; i < tmesh->getHalfEdgeNum(); i++)
//   		{
//   			if(!tmesh->m_pHalfEdge[i].isBoundaryEdge()) 
//   			{
//   				int p1 = tmesh->m_pHalfEdge[i].getVertexIndex(0);
//   				int p2 = tmesh->m_pHalfEdge[i].getVertexIndex(1);
//   				glLineWidth(2.0);
//   				glColor4f(0.0, 0.0, 0.0, 1.0);			//show boundary edge in black
//   				if(tmesh->getVertex_const(p1)->isHole()) 
//   				{
//   					glColor4f(0.0, 0.0, 1.0, 1.0);		//show edge on holes in blue
//   				}
//   				Vector3D v1 = tmesh->getVertex_const(p1)->getPosition();
//   				v1 += shift;
//   				Vector3D v2 = tmesh->getVertex_const(p2)->getPosition();
//   				v2 += shift;
//   				glVertex3d(v1.x, v1.y, v1.z);
//   				glVertex3d(v2.x, v2.y, v2.z);
//   			}
//   		}
//   		glEnd();
//   	}
//  	glEnable(GL_LIGHTING);

	///	draw reference point
	if ( obj == 0 && m_bShowRefPoint && vpMP[0]->pRef >= 0 )
	{
		Vector3D vt = tmesh->getVertex_const(vpMP[0]->pRef)->getPosition();
		if (obj == 0)
			vt = vpMP[0]->posRef;
		vt += shift;
		glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
 		GLUquadric* quadric = gluNewQuadric();
 		gluQuadricDrawStyle(quadric, GLU_FILL);
 		glPushMatrix();
 		glTranslated(vt.x, vt.y, vt.z);
 		gluSphere(quadric, vpMP[obj]->getMesh()->m_edge/4.0, 16, 8);
 		glPopMatrix();
	}

	/// draw handle points
	if (obj == 0)
	{
		for (auto iter = vpMP[0]->mHandles.begin(); iter != vpMP[0]->mHandles.end(); ++iter)
		{
		 	Vector3D vt = iter->second;
		 	vt += shift;
		 	glColor4f(color_handle[0], color_handle[1], color_handle[2], color_handle[3]);
		 	GLUquadric* quadric = gluNewQuadric();
		   	gluQuadricDrawStyle(quadric, GLU_FILL);
		   	glPushMatrix();
		   	glTranslated(vt.x, vt.y, vt.z);
		   	gluSphere(quadric, vpMP[obj]->getMesh()->m_edge/4.0, 16, 8);
		   	glPopMatrix();
		}
	}

	/// draw feature points
	if (obj == 0 && m_bShowFeatures)
	{
		// ---- draw as glPoint ---- //
		glPointSize(10.0);
		glBegin(GL_POINTS);
		for (auto iter = vpMP[0]->vFeatures.begin(); iter != vpMP[0]->vFeatures.end(); ++iter)
		{
			Vector3D vt = tmesh->getVertex_const(iter->index)->getPosition();
			vt += shift;
			glColor4f(featureColors[iter->scale][0], featureColors[iter->scale][1], featureColors[iter->scale][2], featureColors[iter->scale][3]);
			glVertex3d(vt.x, vt.y, vt.z);
		}
		glEnd();
		glPointSize(2.0);

		// ---- draw as gluSphere ---- //
// 		for (auto iter = vpMP[0]->vFeatures.begin(); iter != vpMP[0]->vFeatures.end(); ++iter)
// 		{
// 			Vector3D vt = tmesh->getVertex_const(iter->index)->getPos();
// 			vt += shift;
// 			glColor4f(featureColors[iter->scale][0], featureColors[iter->scale][1], featureColors[iter->scale][2], featureColors[iter->scale][3]);
// 			GLUquadric* quadric = gluNewQuadric();
//   		gluQuadricDrawStyle(quadric, GLU_FILL);
//   		glPushMatrix();
//   		glTranslated(vt.x, vt.y, vt.z);
//   		gluSphere(quadric, vpMP[obj]->mesh->m_edge/4.0, 16, 8);
//   		glPopMatrix();
// 		}
	}

#endif
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}
#endif

void GLMeshWidget::drawMeshExt( const DifferentialMeshProcessor* pMP, const Vector3D& trans, const CQrot& rot, const RenderSettings* renderSettings )
{
	if(!pMP->getMesh()) return;
	const CMesh* tmesh = pMP->getMesh();

	const float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
//  glMateriali(GL_FRONT, GL_SHININESS, 96);

	const Vector3D shift = renderSettings->display_shift;
	const GLfloat *mesh_color = renderSettings->mesh_color;	

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	setupObject(rot, trans);	

	/// draw the mesh
	glPolygonMode(GL_FRONT_AND_BACK, renderSettings->glPolygonMode);

	glPointSize(2.0);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	if (m_bShowSignature && !renderSettings->vDisplaySignature.empty())
	{
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tmesh->getFaceNum(); i++)
		{
			if(!tmesh->m_pFace[i].hasHalfEdge()) continue;
			for (int j = 0; j < 3; j++)
			{
				int pi = tmesh->m_pFace[i].getVertexIndex(j);					 				
				Vector3D norm = tmesh->getVertex_const(pi)->getNormal();
				glNormal3f(norm.x, norm.y, norm.z);	 				
				Vector3D vt = tmesh->getVertex_const(pi)->getPosition();
				vt += shift;	//add some offset to separate object 1 and 2
				glFalseColor(renderSettings->vDisplaySignature[pi], 1.0);	// displaySignature values in [0,1)
				glVertex3f(vt.x, vt.y, vt.z);
			}
		}
		glEnd();
	}
	else // draw mesh in solid color
	{
		glColor4f(mesh_color[0], mesh_color[1], mesh_color[2], mesh_color[3]); 
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tmesh->getFaceNum(); i++)
		{
			if(!tmesh->m_pFace[i].hasHalfEdge()) continue;
			for (int j = 0; j < 3; j++)
			{
				int pi = tmesh->m_pFace[i].getVertexIndex(j);					 				
				Vector3D norm = tmesh->getVertex_const(pi)->getNormal();
				glNormal3f(norm.x, norm.y, norm.z);	 				
				Vector3D vt = tmesh->getVertex_const(pi)->getPosition();
				vt += shift;	//add some offset to separate object 1 and 2
				glVertex3f(vt.x, vt.y, vt.z);
			}
		}
		glEnd();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);

#ifdef DRAW_BOUNDARY
	/// draw boundary edge
	//  	glDisable(GL_LIGHTING);
	//   	if (tmesh->hasBounary())   //highlight boundary edge 
	//   	{
	//   		glBegin(GL_LINES);	
	//   		for(int i = 0; i < tmesh->getHalfEdgeNum(); i++)
	//   		{
	//   			if(!tmesh->m_pHalfEdge[i].isBoundaryEdge()) 
	//   			{
	//   				int p1 = tmesh->m_pHalfEdge[i].getVertexIndex(0);
	//   				int p2 = tmesh->m_pHalfEdge[i].getVertexIndex(1);
	//   				glLineWidth(2.0);
	//   				glColor4f(0.0, 0.0, 0.0, 1.0);			//show boundary edge in black
	//   				if(tmesh->getVertex_const(p1)->isHole()) 
	//   				{
	//   					glColor4f(0.0, 0.0, 1.0, 1.0);		//show edge on holes in blue
	//   				}
	//   				Vector3D v1 = tmesh->getVertex_const(p1)->getPosition();
	//   				v1 += shift;
	//   				Vector3D v2 = tmesh->getVertex_const(p2)->getPosition();
	//   				v2 += shift;
	//   				glVertex3d(v1.x, v1.y, v1.z);
	//   				glVertex3d(v2.x, v2.y, v2.z);
	//   			}
	//   		}
	//   		glEnd();
	//   	}
	//  	glEnable(GL_LIGHTING);
#endif

	//// ----	draw reference point ---- ////
	if ( m_bShowRefPoint && pMP->pRef >= 0 )
	{
		Vector3D vt = tmesh->getVertex_const(pMP->pRef)->getPosition();
		vt = pMP->posRef;
		vt += shift;
		glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
		GLUquadric* quadric = gluNewQuadric();
		gluQuadricDrawStyle(quadric, GLU_FILL);
		glPushMatrix();
		glTranslated(vt.x, vt.y, vt.z);
		gluSphere(quadric, tmesh->m_edge/4.0, 16, 8);
		glPopMatrix();
	}

	//// ---- draw feature points ---- ////
	if (m_bShowFeatures && !pMP->vFeatures.empty())
	{
		///// draw as glPoint
		glPointSize(10.0);
		glBegin(GL_POINTS);
		for (auto iter = pMP->vFeatures.begin(); iter != pMP->vFeatures.end(); ++iter)
		{
			Vector3D vt = tmesh->getVertex_const(iter->index)->getPosition();
			vt += shift;
			int color_index = iter->scale % gFeatureColorNum;
			glColor4f(featureColors[color_index][0], featureColors[color_index][1], featureColors[color_index][2], featureColors[color_index][3]);
			glVertex3d(vt.x, vt.y, vt.z);
		}
		glEnd();
		glPointSize(2.0);

		//// draw as gluSphere 
		// 		for (auto iter = vpMP[0]->vFeatures.begin(); iter != vpMP[0]->vFeatures.end(); ++iter)
		// 		{
		// 			Vector3D vt = tmesh->getVertex_const(iter->index)->getPos();
		// 			vt += shift;
		// 			glColor4f(featureColors[iter->scale][0], featureColors[iter->scale][1], featureColors[iter->scale][2], featureColors[iter->scale][3]);
		// 			GLUquadric* quadric = gluNewQuadric();
		//   		gluQuadricDrawStyle(quadric, GLU_FILL);
		//   		glPushMatrix();
		//   		glTranslated(vt.x, vt.y, vt.z);
		//   		gluSphere(quadric, tmesh->m_edge/4.0, 16, 8);
		//   		glPopMatrix();
		// 		}
	}

	//// ---- draw handle points ---- ////
	if (!pMP->mHandles.empty())
	{
		for (auto iter = pMP->mHandles.begin(); iter != pMP->mHandles.end(); ++iter)
		{
			Vector3D vt = iter->second;
			vt += shift;
			glColor4f(color_handle[0], color_handle[1], color_handle[2], color_handle[3]);
			GLUquadric* quadric = gluNewQuadric();
			gluQuadricDrawStyle(quadric, GLU_FILL);
			glPushMatrix();
			glTranslated(vt.x, vt.y, vt.z);
			gluSphere(quadric, tmesh->m_edge/4.0, 16, 8);
			glPopMatrix();
		}
	}
	
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void GLMeshWidget::drawLegend(QPainter* painter)
{
	painter->setRenderHint(QPainter::Antialiasing);
	
	int xBegin = width()/2 - 128;

	for (int i = 0; i < 255; i++)
	{
		QColor col(255*falseColorMap.RedMap[i], 255*falseColorMap.GreenMap[i], 255*falseColorMap.BlueMap[i], 255);
//		painter->setBrush(QBrush(col, Qt::SolidPattern));
//		painter->drawRect(xBegin + i*2, height()-50, 2, 25);
		painter->setPen(QPen(col, 1, Qt::SolidLine));
		painter->drawLine(QPointF(xBegin+i, height()-50), QPointF(xBegin+i, height()-25));
	}
	painter->setPen(QPen(Qt::black, Qt::SolidLine));
	painter->drawText(xBegin, height() - 70, 128, 12, Qt::AlignLeft, QString::number(vpMP[0]->sigMin));
	painter->drawText(xBegin + 128, height()-70, 128, 12, Qt::AlignRight, QString::number(vpMP[0]->sigMax));
}

// void GLMeshWidget::showEvent( QShowEvent *event )
// {
// 	Q_UNUSED(event);
// }

void GLMeshWidget::setupViewport( int width, int height )
{
	GLdouble ar = GLdouble(width) / height;	//aspect ratio

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

//	glFrustum(-ar, ar, -1.0, 1.0, 4.0, 15.0);
// 	GLdouble clipX = g_myNear * tan(g_myAngle/2.0/180.0 * PI), 
// 		     clipY = clipX / ar;
// 	glFrustum(-clipX, -clipY, clipX, clipY, g_myNear, g_myFar);

	gluPerspective(g_myAngle, ar, g_myNear, g_myFar);
//	glMatrixMode(GL_MODELVIEW);
}

void GLMeshWidget::paintEvent( QPaintEvent *event )
{
	// To achieve 2D graphics and 3d OpenGL overlay, we have to implement paintEvent instead of relying on paintGL()

	QPainter painter(this);
	drawGL();

	if (m_bShowLegend && vSettings[0].vDisplaySignature.empty())
 	drawLegend(&painter);
}

bool GLMeshWidget::glPick( int x, int y, Vector3D& _p )
{
	GLdouble  modelview[16], projection[16];
	GLint     viewport[4];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
	CQrot& rot = ObjRot1;
	Vector3D& trans = ObjTrans1;
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

	if (z != 1.0f)
	{
		_p = Vector3D(pos[0], pos[1], pos[2]);
		return true;
	}

	return false;

}


