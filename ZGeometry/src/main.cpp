#include "QZGeometry.h"
#include <QtGui/QApplication>
#include <QtGui/QLabel>
#include <QGLFormat>
#include <SimpleConfigLoader.h>
#include "OutputHelper.h"

SimpleConfigLoader g_configMgr("configs.txt");
OutputHelper qout;
QString qformat;

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	QGLFormat fmt;
	fmt.setProfile(QGLFormat::CompatibilityProfile);
//	fmt.setVersion(4, 3);
	QGLFormat::setDefaultFormat(fmt);

	QZGeometryWindow w;
	if (!w.initialize()) exit(-1);
	w.show();
	return a.exec();
}
