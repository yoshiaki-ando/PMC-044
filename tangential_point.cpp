#include <cmath>
#include <Vector3d.h>
#include "pmc_simulation.h"

/*
 * yz平面において、角度αだけ(+xに対して右ねじ方向に)回転
 *
 * [x']   [1   0      0  ][x]
 * [y'] = [0 cosα -sinα][y]
 * [z']   [0 sinα  cosα][z]
 */

AndoLab::Vector3d <double> rotate_alpha(
    AndoLab::Vector3d <double> r, const double alpha){

  return AndoLab::Vector3d <double>
  ( r.x(),
      std::cos(alpha)*r.y() - std::sin(alpha)*r.z(),
      std::sin(alpha)*r.y() + std::cos(alpha)*r.z() );
}

/*
 * zx平面において、角度θだけ(+yに対して右ねじ方向に)回転
 *
 * [x']   [ cosθ 0 sinθ][x]
 * [y'] = [  0    1   0  ][y]
 * [z']   [-sinθ 0 cosθ][z]
 */

AndoLab::Vector3d <double> rotate_theta(
    AndoLab::Vector3d <double> r, const double theta){

  return AndoLab::Vector3d <double>
  ( std::cos(theta)*r.x() + std::sin(theta)*r.z(),
      r.y(),
      -std::sin(theta)*r.x() + std::cos(theta)*r.z() );
}

AndoLab::Vector3d <double> tangential_point(double alt, double alpha){
  double R = Radius_of_Earth + alt;
  double cos_b = R / Rgeo;
  double sin_b = std::sqrt( 1.0 - cos_b*cos_b );
  AndoLab::Vector3d <double> r(R*cos_b, 0.0, R*sin_b );

  return rotate_alpha( r, alpha );
//  return AndoLab::Vector3d <double> (
//      R*cos_b, -R*sin_b*sin(alpha), R*sin_b*cos(alpha) );
}
