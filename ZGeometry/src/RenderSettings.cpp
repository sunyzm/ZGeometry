#include "RenderSettings.h"
#include "global.h"

RenderSettings::RenderSettings() :
    displayType(Mesh), selected(false),
    glPolygonMode(0x1B02), display_shift(0, 0, 0),
    obj_rot(1, 0, 0, 0), obj_trans(0, 0, 0),
    mActiveColorSignatureName(StrAttrColorPreset)
{

}
