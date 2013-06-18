#include "QZGeometry.h"
#include <QtGui/QApplication>
#include <QtGui/QLabel>
#include <QGLFormat>
#include <SimpleConfigLoader.h>
#include "OutputHelper.h"
#include "common.h"

#if _MSC_VER < 1600
#error VC++ 2010 or later required.
#endif

SimpleConfigLoader g_configMgr("configs.txt");
OutputHelper qout;
GeometryTask g_task = TASK_REGISTRATION;

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

// 	QGLFormat fmt;
// 	fmt.setProfile(QGLFormat::CompatibilityProfile);
// 	QGLFormat::setDefaultFormat(fmt);
	
	g_task = (GeometryTask)g_configMgr.getConfigValueInt("GEOMETRY_TASK");

	QZGeometryWindow w;
	if (!w.initialize()) exit(-1);
	w.show();
	return a.exec();
}
