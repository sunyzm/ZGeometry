#pragma once

#define DEFAULT_EIGEN_SIZE 500
#define DEFAULT_DEFORM_RING 5 
#define LOAD_MHB_CACHE 1
#define MIN_HK_TIMESCALE 1e-2
#define DEFUALT_HK_TIMESCALE 40.0
#define MAX_HK_TIMESCALE 2000.0
#define PARAMETER_SLIDER_CENTER 50

#define DR_THRESH_INCREMENT 0.00001
#define MATCHING_THRESHOLD 0.002

enum GeometryTask {TASK_VIEWING, TASK_REGISTRATION, TASK_EDITING};

extern const char* g_meshListName;
extern const int g_preload_meshes;