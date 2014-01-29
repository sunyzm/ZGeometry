#include "global.h"

OutputHelper qout;
SimpleConfigLoader g_configMgr("config.txt");
GeometryTask g_task = TASK_REGISTRATION;
ZGeom::MatlabEngineWrapper g_engineWrapper(256);	
