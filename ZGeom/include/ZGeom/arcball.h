#ifndef ZGEOM_ARCBALL_H
#define ZGEOM_ARCBALL_H
#include "quat.h"
#include "Vec2.h"
#include "Vec3.h"

class CArcball
{
public:
  CArcball(){};

  CArcball( int width, int height, int ox, int oy );
  CQrot update( int nx, int ny );

private:
  void _plane2sphere( const ZGeom::Vec2d & v, ZGeom::Vec3d & r );

  ZGeom::Vec3d m_position;
  double   m_radius;
  ZGeom::Vec2d m_center;
};

#endif