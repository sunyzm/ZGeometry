#ifndef GLMESHWIDGET_H
#define GLMESHWIDGET_H

#include <QGLWidget>

class GLMeshWidget : public QGLWidget
{
	Q_OBJECT

public:
	GLMeshWidget(QWidget *parent = 0);
	~GLMeshWidget();
protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();

private:
	void draw();

	GLfloat rotationX, rotationY, rotationZ;
	QColor faceColors[4];
};

#endif // GLMESHWIDGET_H
