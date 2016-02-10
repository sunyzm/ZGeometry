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

void loadConfig()
{
    /* read in configuration parameters from g_configMgr */
    g_configMgr.getConfigValueInt("GEOMETRY_TASK", (int&)g_task);
    g_configMgr.getConfigValueInt("LOAD_MHB_CACHE", gSettings.LOAD_MHB_CACHE);
    g_configMgr.getConfigValueDouble("PARAMETER_SLIDER_CENTER", gSettings.PARAMETER_SLIDER_CENTER);
    g_configMgr.getConfigValueDouble("DEFUALT_HK_TIMESCALE", gSettings.DEFAULT_HK_TIMESCALE);
    g_configMgr.getConfigValueInt("DEFAULT_EIGEN_SIZE", gSettings.DEFAULT_EIGEN_SIZE);
    g_configMgr.getConfigValueInt("SCALE_TO_UNIT", gSettings.INPUT_SCALE_TO_UNIT);
    g_configMgr.getConfigValueInt("PARITION_SIZE", gSettings.MESH_PARTITION_SIZE);

}

int main(int argc, char *argv[])
{
    loadConfig();
	
    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setProfile(QSurfaceFormat::CompatibilityProfile);
    QSurfaceFormat::setDefaultFormat(format);

    QApplication a(argc, argv);
    QZGeometryWindow w;
    QObject::connect(&w, SIGNAL(displayQtVersion()), &a, SLOT(aboutQt()));

    CStopWatch timer;
    timer.startTimer();
    g_engineWrapper.open();
    timer.stopTimer("-- Matlab Initialization time: ", " --");

    if (!w.initialize(mesh_list_name)) std::exit(-1);

    w.show();
	return a.exec();
}
