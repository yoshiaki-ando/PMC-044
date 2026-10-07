/*
 * geocoordinate.h
 *
 * 緯度、経度、高度
 * 「ひまわり座標(東経140°をφ=0とした、[m]の単位を持つ座標)」との変換
 *
 *  Created on: 2019/07/10
 *      Author: ando
 */

#ifndef GEOCOORDINATE_H_
#define GEOCOORDINATE_H_

#include <cmath>
#include <Vector3d.h>

constexpr double Meter2killometer { 1.0e-3 };

class Geocoordinate{
private:
  double pLatitude;  /* [deg] */
  double pLongitude; /* [deg] */
  double pAltitude;  /* [m] */
  AndoLab::Vector3d <double> pR;

  AndoLab::Vector3d <double> projection_on_yz(void);

public:
  Geocoordinate(void):
    pLatitude(0.0), pLongitude(0.0), pAltitude(0.0) { }

  Geocoordinate(AndoLab::Vector3d <double> r);
  Geocoordinate(const double Latitude, const double Longitude, const double Altitude);

  /* アクセサ */
  void set(AndoLab::Vector3d <double> r);
  void set(const double Latitude, const double Longitude, const double Altitude);
  double latitude(void){ return pLatitude; }
  double longitude(void){ return pLongitude; }
  double east_longitude(void); /* 全て東経で述べた経度 */
  double altitude(void){ return pAltitude; }
  double altitude_in_km(void){ return pAltitude*Meter2killometer; }
  AndoLab::Vector3d <double> r(void){ return pR; }

  /* 衛星からみた角度αに直す */
  double alpha(void);
};



#endif /* GEOCOORDINATE_H_ */
