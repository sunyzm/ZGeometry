#include <cstdlib>
#include <QApplication>
#include <ppl.h>
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
	
    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setProfile(QSurfaceFormat::OpenGLContextProfile::CompatibilityProfile);
    QSurfaceFormat::setDefaultFormat(format);

    QApplication a(argc, argv);
    QZGeometryWindow w;
    QObject::connect(&w, SIGNAL(displayQtVersion()), &a, SLOT(aboutQt()));

    CStopWatch timer;
    timer.startTimer();
    g_engineWrapper.open();
    if (!w.initialize(mesh_list_name)) std::exit(-1);
    timer.stopTimer("-- Initialization time: ", " --");
    
    w.show();
	return a.exec();
}
