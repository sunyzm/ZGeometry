#ifndef ZGEOM_COLOR_MAP_H
#define ZGEOM_COLOR_MAP_H
#include "Color.h"

namespace ZGeom {

enum ColorMapType { 
    CM_PSEUDO = 0, CM_PARULA , CM_JET, CM_COOL, CM_HOT, COLOR_MAP_COUNT
};
    
class ColorMap {
public:
    static Colorf falseColor(float gray, float alpha = 1.0, ColorMapType cmt = CM_PSEUDO);
    static double *getColorMap(ColorMapType map_type);
    static double parula[];
    static double jet[];
};

}   // end of namespace

#endif