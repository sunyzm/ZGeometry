#include "glmeshwidget.h"
#include <util/OutputHelper.h>

static OutputHelper qout;

GLMeshWidget::GLMeshWidget(QWidget *parent)
	: QGLWidget(parent)
{
	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));
	faceColors[0] = Qt::red;
	faceColors[1] = Qt::green;
	faceColors[2] = Qt::blue;
	faceColors[3] = Qt::yellow;

	rotationX = -21.0;
	rotationY = -57.0;
	rotationZ = 0.0;
}

GLMeshWidget::~GLMeshWidget()
{

}

void GLMeshWidget::initializeGL()
{
//	qglClearColor(Qt::black);
	glClearColor(.0f, .0f, .0f, 1.0f);
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
}

void GLMeshWidget::resizeGL( int width, int height )
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLfloat x = GLfloat(width) / height;
	glFrustum(-x, x, -1.0, 1.0, 4.0, 15.0);
	glMatrixMode(GL_MODELVIEW);
}

void GLMeshWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw();
}

void GLMeshWidget::draw()
{
	static const GLfloat P1[3] = {0, -1, 2};
	static const GLfloat P2[3] = {1.732, -1, -1};
	static const GLfloat P3[3] = {-1.732, -1, -1};
	static const GLfloat P4[3] = {0, 2, 0};

	static const GLfloat *const coords[4][3] = {{P1, P2, P3}, {P1,P3,P4}, {P1,P4,P2}, {P2,P4,P3}};

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0,0,-10);
	glRotatef(rotationX, 1, 0, 0);
	glRotatef(rotationY, 0, 1, 0);
	glRotatef(rotationZ, 0, 0, 1);

	for (int i = 0; i < 4; ++i)
	{
		glLoadName(i);
		glBegin(GL_TRIANGLES);
		qglColor(faceColors[i]);
		for (int j = 0; j < 3; ++j)
		{
			glVertex3f(coords[i][j][0], coords[i][j][1], coords[i][j][2]);
		}
		glEnd();
	}

	qout.output("Welcome from GLMesh!", OUT_CONSOLE);
}
