#include <cstdlib>
#include <QApplication>
#include "QZGeometry.h"
#include "global.h"

#ifdef NDEBUG
std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME");
#else NDEBUG
std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME_DEBUG");
#endif

int main(int argc, char *argv[])
{
	g_configMgr.getConfigValueInt("GEOMETRY_TASK", (int&)g_task);
	CStopWatch timer;
	timer.startTimer();
	g_engineWrapper.open();     
	timer.stopTimer("-- Matlab engine opening time: ", " --");

	QApplication a(argc, argv);
	QZGeometryWindow w;
	if (!w.initialize(mesh_list_name)) std::exit(-1);

	w.show();
	return a.exec();
}
