#include "glmeshwidget.h"
#include <util/OutputHelper.h>
#include <util/util.h>
#include <gl/glu.h>

const GLfloat color1[4] = {0.53, 0.70, 0.93, 1.0};
const GLfloat color2[4] = {0.99, 0.73, 0.62, 1.0}; //{0.63,0.78,0.63,1.0};

extern OutputHelper qout;
extern QString qformat;
//extern int g_objSelect;
Qt::MouseButton gButton;


void glFalseColor(float v, float p)
{
	int floor = v * 255.0;
	glColor4f(FalseColorMap::RedMap[floor], FalseColorMap::GreenMap[floor], FalseColorMap::BlueMap[floor], p);
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

	bShowRefPoint = false;
	FalseColorMap::BuildLUT();

	vSettings.resize(2, DisplaySettings());
}

GLMeshWidget::~GLMeshWidget()
{

}

void GLMeshWidget::addMesh(MeshProcessor* pmp)
{
	vpMP.push_back(pmp);
}

void GLMeshWidget::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton)
	{
		gButton = Qt::LeftButton;
		g_arcball = CArcball(width(), height(), event->x() - width()/2, height()/2 - event->y());
	}
	else if (event->button() == Qt::MidButton)
	{
		gButton = Qt::MidButton;
		g_startx = event->x();
		g_starty = event->y();
	}
	else if (event->button() == Qt::RightButton)
	{
		gButton = Qt::RightButton;
		g_startx = event->x();
		g_starty = event->y();
	}
}

void GLMeshWidget::mouseMoveEvent(QMouseEvent *event)
{
	Vector3D trans;
	CQrot    rot;

	int win_width = this->width(), win_height = this->height();
	int x = event->x(), y = event->y();

	if (gButton == Qt::LeftButton ) 
	{
		rot = g_arcball.update( x - win_width/2, win_height - y - win_height/2);
		if (vSettings[0].selected) 
			ObjRot1 = rot * ObjRot1;
		if (vSettings[1].selected) 
			ObjRot2 = rot * ObjRot2;
	}
	else if (gButton == Qt::MidButton) 
	{
		float scale = 3.0 * vpMP[0]->mesh->m_bBox.x / win_height;
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
		float scale = 5.0 * vpMP[0]->mesh->m_bBox.y / win_height;
		trans =  Vector3D(0, 0, scale * (g_starty - y));
		g_startx = x;
		g_starty = y;
		if (vSettings[0].selected)
			ObjTrans1 = ObjTrans1 + trans;
		if (vSettings[1].selected) 
			ObjTrans2 = ObjTrans2 + trans;
	}

	updateGL();
}

void GLMeshWidget::wheelEvent(QWheelEvent *event)
{
	if (event->modifiers() & Qt::ControlModifier)
	{
		int numSteps = event->delta();
		float scale = 3.0 * vpMP[0]->mesh->m_bBox.x / this->height();
		Vector3D trans =  Vector3D(0, 0, scale * numSteps);
		
		if (vSettings[0].selected) 
			ObjTrans1 = ObjTrans1 + trans;
		if (vSettings[1].selected)
			ObjTrans2 = ObjTrans2 + trans;		
	}

	updateGL();
}

void GLMeshWidget::initializeGL()
{
	GLfloat diffuse[] = {1, 1, 1, 1};
	GLfloat global_ambient[] = {.2, .2, .2, 1};
	GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat position[] = {.0,  .0, 1, 0.0};

	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);      
	glEnable(GL_DEPTH_TEST);
	glClearColor(1,1,1,0);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable (GL_POLYGON_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

	glPolygonMode(GL_FRONT, GL_FILL);
	glEnable (GL_BLEND); 
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
	//glDisable(GL_DEPTH_TEST);

	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);

	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glLightfv(GL_LIGHT1, GL_POSITION, position);
}

void GLMeshWidget::fieldView( const Vector3D &center, const Vector3D &bbox )
{
	g_EyeZ = 15.0 * (float)bbox.y;
	g_myFar = 100.0 * (float)bbox.y;
	g_myNear = 0.01 * g_myFar;
	float len = (bbox.y > bbox.x) ? (float)bbox.y : (float)bbox.x;
	g_myAngle = 2.0 * atan2(len, g_EyeZ);
	g_myAngle = (g_myAngle * 180.0) / PI + 5.0;

	qout.output(qformat.sprintf("(%f,%f,%f) - (%f,%f,%f)", center.x, center.y, center.z, bbox.x, bbox.y, bbox.z));
}

void GLMeshWidget::resizeGL( int width, int height )
{
	GLdouble ar = GLdouble(width) / height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

//	glFrustum(-ar, ar, -1.0, 1.0, 4.0, 15.0);

// 	GLdouble clipX = g_myNear * tan(g_myAngle/2.0/180.0 * PI), 
// 		     clipY = clipX / ar;
// 	glFrustum(-clipX, -clipY, clipX, clipY, g_myNear, g_myFar);
	gluPerspective(g_myAngle, ar, g_myNear, g_myFar);
	glMatrixMode(GL_MODELVIEW);
}

void GLMeshWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw();
}

void GLMeshWidget::draw()
{
	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
	
// 	if (vSettings[0].displayType == DisplaySettings::Mesh)
// 		drawMesh(vpMP[0]->mesh, ObjRot1, ObjTrans1, color1);
	
	if (vSettings[0].displayType == DisplaySettings::Mesh || vSettings[0].displayType == DisplaySettings::Signature)
	{
		drawMeshExt(0);
		drawMeshExt(1);
	}
}

void GLMeshWidget::setupObject(const CQrot& qrot, const Vector3D& trans)
{
	glTranslated(trans.x, trans.y, trans.z);
	double rot[16];
	qrot.convert( rot );
	glMultMatrixd(( GLdouble*)  rot );
}

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
			if(!tmesh->m_pFace[i].m_piEdge) continue;
			for (int j = 0; j < 3; j++)
			{
				glColor4f(color[0], color[1], color[2], color[3]); 
				int pi = tmesh->m_pFace[i].m_piVertex[j];
				Vector3D norm = tmesh->m_pVertex[pi].getNormal();
				glNormal3f(norm.x, norm.y, norm.z);
				Vector3D vt = tmesh->m_pVertex[pi].m_vPosition;
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
				Vector3D v1 = tmesh->m_pVertex[p1].m_vPosition;
				//v1 -= tmesh->m_Center;
				Vector3D v2 = tmesh->m_pVertex[p2].m_vPosition;
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

void GLMeshWidget::drawMeshExt( int obj )
{
	if (obj >= vpMP.size() || obj < 0) return;	
	if(!vpMP[obj]->mesh) return;

	CMesh* tmesh = vpMP[obj]->mesh;
	CQrot rot = (obj == 0) ? ObjRot1 : ObjRot2;
	Vector3D trans = (obj == 0) ? ObjTrans1 : ObjTrans2;
	const GLfloat *color = (obj == 0) ? color1 : color2;	

	float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
//  glMateriali(GL_FRONT, GL_SHININESS, 96);

	Vector3D shift = Vector3D(tmesh->m_bBox.x/2, 0, 0);
	if (obj > 0) shift.x = -shift.x;

	glPushMatrix();
	setupObject(rot, trans);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);
	
	if (vSettings[obj].displayType == DisplaySettings::Mesh || vpMP[obj]->vDisplaySignature.empty())
	{	// just display mesh in single color
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tmesh->getFaceNum(); i++)
		{
			if(!tmesh->m_pFace[i].m_piEdge) continue;
			for (int j = 0; j < 3; j++)
			{
				int pi = tmesh->m_pFace[i].m_piVertex[j];
				glColor4f(color[0], color[1], color[2], color[3]); 
				Vector3D norm = tmesh->m_pVertex[pi].getNormal();
				glNormal3f(norm.x, norm.y, norm.z);
				Vector3D vt = tmesh->m_pVertex[pi].m_vPosition;
				vt -= shift;
				glVertex3f(vt.x, vt.y, vt.z);
			}
		}
		glEnd();
	}
	else 
	{	
		//display signature value in false color
	 	glBegin(GL_TRIANGLES);
	 	for (int i = 0; i < tmesh->getFaceNum(); i++)
	 	{
	 		if(!tmesh->m_pFace[i].m_piEdge) continue;
	 		for (int j = 0; j < 3; j++)
	 		{
	 			int pi = tmesh->m_pFace[i].m_piVertex[j];
	 			float scaleVal = vpMP[obj]->vDisplaySignature[pi];
	 			glFalseColor(scaleVal, 1.0f);	 				
	 			Vector3D norm = tmesh->m_pVertex[pi].getNormal();
	 			glNormal3f(norm.x, norm.y, norm.z);	 				
	 			Vector3D vt = tmesh->m_pVertex[pi].m_vPosition;
	 			vt -= shift;
	 			glVertex3f(vt.x, vt.y, vt.z);
	 		}
	 	}
	 	glEnd();
	}
	
	glDisable(GL_POLYGON_OFFSET_FILL);
	
	glDisable(GL_LIGHTING);
	if (tmesh->hasBounary())   //highlight boundary edge 
	{

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
				Vector3D v1 = tmesh->m_pVertex[p1].m_vPosition;
				v1 -= shift;
				Vector3D v2 = tmesh->m_pVertex[p2].m_vPosition;
				v2 -= shift;
				glVertex3d(v1.x, v1.y, v1.z);
				glVertex3d(v2.x, v2.y, v2.z);
			}
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);

//	control display of reference point
	if (vpMP[0]->pRef >= 0 && bShowRefPoint && obj == 0)
	{
		Vector3D vt = tmesh->m_pVertex[vpMP[0]->pRef].m_vPosition;
		if (obj == 0)
			vt = vpMP[0]->posRef;
		vt -= shift;
		glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
		GLUquadric* quadric = gluNewQuadric();
		gluQuadricDrawStyle(quadric, GLU_FILL);
		glPushMatrix();
		glTranslated(vt.x, vt.y, vt.z);
		gluSphere(quadric, vpMP[obj]->mesh->m_edge, 8, 8);
		glPopMatrix();
	}

	glPopMatrix();

}
