#include "glmeshwidget.h"
#include <util/OutputHelper.h>
#include <util/util.h>
#include <gl/glu.h>

extern OutputHelper qout;
const GLfloat color1[4] = {0.53, 0.70, 0.93, 1.0};
const GLfloat color2[4] = {0.99, 0.73, 0.62, 1.0}; //{0.63,0.78,0.63,1.0};
QString qformat;

GLMeshWidget::GLMeshWidget(QWidget *parent)
	: QGLWidget(parent)
{
	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer | QGL::Rgba));

	g_EyeZ = 10.0;
	ObjRot1 = ObjRot2 = CQrot(1,0,0,0);
	ObjTrans1 = ObjTrans2 = Vector3D(0,0,0);
	g_myNear = 1.0;
	g_myFar = 100.0;
	g_myAngle = 40.0;
}

GLMeshWidget::~GLMeshWidget()
{

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
//	glClearColor(1,1,1,0);
	draw();
}

void GLMeshWidget::draw()
{
	glLoadIdentity();
	gluLookAt(0, 0, g_EyeZ, 0, 0, 0, 0, 1, 0);
	drawMesh(mp[0].mesh, ObjRot1, ObjTrans1, color1);
}

void GLMeshWidget::setMesh( CMesh* cm, int i/* = 0*/ )
{
	mp[i].mesh = cm;
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
		for (int i = 0; i < tmesh->m_nFace; i++)
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

//	glEnable(GL_LIGHTING);
	glPopMatrix();
}