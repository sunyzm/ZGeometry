#include <cstdlib>
#include <QApplication>
#include "QZGeometry.h"
#include "global.h"

#ifdef NDEBUG
std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME");
#else
std::string mesh_list_name = g_configMgr.getConfigValue("MESH_LIST_NAME_DEBUG");
#endif

int main(int argc, char *argv[])
{
	/* read in configuration parameters from g_configMgr */
	g_configMgr.getConfigValueInt("GEOMETRY_TASK", (int&)g_task);
	g_configMgr.getConfigValueInt("LOAD_MHB_CACHE", gSettings.LOAD_MHB_CACHE);
	g_configMgr.getConfigValueDouble("PARAMETER_SLIDER_CENTER", gSettings.PARAMETER_SLIDER_CENTER);
	g_configMgr.getConfigValueDouble("DEFUALT_HK_TIMESCALE", gSettings.DEFUALT_HK_TIMESCALE);
	g_configMgr.getConfigValueInt("DEFAULT_EIGEN_SIZE", gSettings.DEFAULT_EIGEN_SIZE);
	
    CStopWatch timer;
    timer.startTimer();
    g_engineWrapper.open();
    timer.stopTimer("-- Matlab engine opening time: ", " --");

	QApplication a(argc, argv);
	QZGeometryWindow w;
	QObject::connect(&w, SIGNAL(displayQtVersion()), &a, SLOT(aboutQt()));
	if (!w.initialize(mesh_list_name)) std::exit(-1);
	w.show();
	return a.exec();
}
