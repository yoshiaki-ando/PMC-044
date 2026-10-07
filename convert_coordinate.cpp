/*
 * convert_coordinate.cpp
 *
 *  Created on: 2024/03/04
 *      Author: ando
 *
 * 緯度・経度から、計算に使う座標(ひまわり方向を +x方向)への変換
 *
 * 緯度・経度・高度 0m の位置を r0 とする。
 */


#include <Vector3d.h>
#include "pmc_simulation.h"

void convert_coordinate(const double latitude, const double longitude,
    AndoLab::Vector3d <double> &r0){

  /* ひまわり140.7°をφ=0とした座標系での phi */
  double geo_phi { phi( lat2theta(latitude) ) };

  /* ひまわりの東経より低い場合はφ<0 */
  if ( (longitude < Longitude_of_himawari) &&
      (longitude > 180.0 - Longitude_of_himawari) ){
    geo_phi *= -1.0;
  }

  double theta = lat2theta( latitude );

  r0.set( Radius_of_Earth, theta, geo_phi, AndoLab::coordinate::Spherical);
}

