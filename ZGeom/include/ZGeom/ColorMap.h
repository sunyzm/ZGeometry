#ifndef ZGEOM_COLOR_MAP_H
#define ZGEOM_COLOR_MAP_H

namespace ZGeom {
    enum ColorMapType { CM_PARULA = 0, CM_JET, CM_COOL, CM_HOT, CM_COUNT };
    
    class ColorMap {
    public:
        static double parula[];
        static double jet[];
    };
}
#endif