#include "global.h"

OutputHelper qout;
SimpleConfigLoader g_configMgr("config.txt");
GeometryTask g_task = TASK_REGISTRATION;
ZGeom::MatlabEngineWrapper g_engineWrapper(256);	
GlobalSettings gSettings;


double inpainting_error_curving_max = 0.03;

GlobalSettings::GlobalSettings()
{
	DEFAULT_EIGEN_SIZE      = -1;
	DEFAULT_DEFORM_RING     = 5;
	LOAD_MHB_CACHE          = 0;
	MIN_HK_TIMESCALE        = 1;
	DEFUALT_HK_TIMESCALE    = 40.0;
	MAX_HK_TIMESCALE        = 2000.0;
	PARAMETER_SLIDER_CENTER = 50;
	DR_THRESH_INCREMENT     = 0.00001;
}
