#include <cstdlib>
#include <QtGui/QApplication>
#include "QZGeometry.h"
#include "global.h"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	
	g_task = static_cast<GeometryTask>(g_configMgr.getConfigValueInt("GEOMETRY_TASK"));

	QZGeometryWindow w;

	if (!w.initialize()) std::exit(-1);

	w.show();
	return a.exec();
}
