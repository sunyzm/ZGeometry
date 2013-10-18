#include <algorithm>
#include <cassert>
#include <tuple>
#include "RenderSettings.h"
#include "global.h"

using namespace std;

RenderSettings::RenderSettings() : mesh_color(preset_mesh_colors[0]), displayType(Mesh), 
	showFeatures(false), showRefPoint(false), 
	showColorSignature(false), selected(false), 
	glPolygonMode(Z_FILL), display_shift(0, 0, 0),
	obj_rot(1,0,0,0), obj_trans(0,0,0)
{
	mColorSignatureName = StrColorPreset;
}
