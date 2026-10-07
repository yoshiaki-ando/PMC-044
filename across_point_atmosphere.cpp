#include <cmath>
#include <array>
#include <Vector3d.h>
#include "Fixed_parameter.h"
#include "pmc_simulation.h"

/* ひまわりから r の点へ向かう直線で大気圏頂上と交差する点を求める
 *
 * return_vector[0]が近い点、[1]が遠い点
 */
std::array<AndoLab::Vector3d <double>, 2> Across_point_atmosphere(AndoLab::Vector3d <double> r){
  AndoLab::Vector3d <double> r_geo(Rgeo, 0.0, 0.0); /* ひまわりの位置 */
//  AndoLab::Vector3d <double> re = r - r_geo; /* ひまわりから、その点に向かうベクトル */
//  AndoLab::Vector3d <double> ren = re.n(); /* その単位ベクトル */

  /* ひまわりから r の点へ向かう直線で大気圏頂上と交差する点
   * r_geo + t ren = R0 r(θ、φ) となる t */
//  double t[2];
//  t[0] = -Rgeo*ren.x()
//      - std::sqrt( Rgeo*Rgeo * (ren.x()*ren.x() - 1.0) + R_Atmosphere*R_Atmosphere );
//  t[1] = -Rgeo*ren.x()
//      + std::sqrt( Rgeo*Rgeo * (ren.x()*ren.x() - 1.0) + R_Atmosphere*R_Atmosphere );
//  AndoLab::Vector3d <double> r1 = 1.0e-3*(r_geo+t1*ren);
//  AndoLab::Vector3d <double> r2 = 1.0e-3*(r_geo+t2*ren);
  std::array<AndoLab::Vector3d <double>, 2> r_atmos;
//  for(int i = 0; i < 2; i++){
//    r_atmos[i] = r_geo + t[i]*ren;
//  }

  AndoLab::Vector3d <double> re = r_geo - r; /* その点から、ひまわりに向かうベクトル */
  r_atmos[0] = cross_point_at_altitude(r, re.n(), Altitude_of_Atmosphere);
  r_atmos[1] = cross_point_at_altitude(r, -1.0*re.n(), Altitude_of_Atmosphere);

  return r_atmos;
}
