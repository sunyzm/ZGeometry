#include "arcball.h"

using ZGeom::Vec2d;
using ZGeom::Vec3d;

CArcball::CArcball( int win_width, int win_height, int ox, int oy )
{
	m_radius = (win_width < win_height ) ? win_width/2 : win_height/2;
	m_center = ZGeom::Vec2d( win_width/2, win_height/2 );

	ZGeom::Vec2d p(ox, oy);
	_plane2sphere(p, m_position);
}

void CArcball::_plane2sphere(const ZGeom::Vec2d & p, Vec3d & q)
{

  ZGeom::Vec2d f = p;
  f /= m_radius;
  double l = f.length();

  if( l > 1.0 ) {
	  q = Vec3d( f.x/l, f.y/l, 0);
	  return;
  }

  double fz = sqrt( 1.0 - l*l );
  q = Vec3d(f.x, f.y, fz);
}

CQrot CArcball::update( int nx, int ny )
{
	Vec3d position;
	_plane2sphere(ZGeom::Vec2d(nx, ny), position);
	Vec3d cp = ZGeom::cross<double>(m_position, position);
	CQrot r(m_position.dot(position), cp.x, cp.y, cp.z);
	m_position = position;

	return r;
}