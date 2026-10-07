/*
 * cross_point_at_altitude.cpp
 *
 *  Created on: 2022/12/21
 *      Author: ando
 */
#include <cmath>
#include <Vector3d.h>
#include "pmc_simulation.h"

/* rの点から、d方向へ進んで高度 H (原点からの距離 = 地球の半径 + H)
 * との交点を求める
 */
AndoLab::Vector3d <double> cross_point_at_altitude(
    AndoLab::Vector3d <double> r,
    AndoLab::Vector3d <double> d, /* 単位ベクトルとすること */
    const double H /* [m] 高度 */
    ){

  const double R_H = Radius_of_Earth + H; /* 原点からの距離 */
  double r_dot_d = r%d;
  double r_abs = r.abs();
  double k = -1.0*r_dot_d +
      std::sqrt( r_dot_d*r_dot_d - r_abs*r_abs + R_H*R_H );
  return r + k*d;
}



