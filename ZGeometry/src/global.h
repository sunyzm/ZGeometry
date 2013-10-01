#pragma once

#include <ZUtil/SimpleConfigLoader.h>
#include "OutputHelper.h"

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};
enum DeformType {Simple, Gradient, Laplace, BiLaplace, Shell, SGW};

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;
extern GeometryTask g_task;