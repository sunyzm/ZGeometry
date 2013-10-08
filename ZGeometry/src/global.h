#pragma once
#include <string>
#include <ZUtil/SimpleConfigLoader.h>
#include "OutputHelper.h"

enum GeometryTask {TASK_VIEWING = 0, TASK_REGISTRATION = 1, TASK_EDITING = 2};
enum DeformType {Simple, Gradient, Laplace, BiLaplace, Shell, SGW};

const std::string StrColorGaussCurvature = "vert_color_gauss_curvature";
const std::string StrColorMeanCurvature  = "vert_color_mean_curvature";
const std::string StrColorEigenFunction  = "vert_color_eig_function";
const std::string StrColorHKS			 = "vert_color_hks";
const std::string StrColorHK			 = "vert_color_hk";
const std::string StrColorHeat			 = "vert_color_heat";
const std::string StrColorMHWS			 = "vert_color_mhws";
const std::string StrColorPreset		 = "vert_color_preset";

extern OutputHelper qout;
extern SimpleConfigLoader g_configMgr;
extern GeometryTask g_task;