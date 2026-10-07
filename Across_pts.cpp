#include <iostream>
#include <Vector3d.h>

#include "Fixed_parameter.h"
#include "Across_pts.h"

Across_pts::Across_pts(AndoLab::Vector3d <double> r2, AndoLab::Vector3d <double> r1,
    double Altitude){

  AndoLab::Vector3d <double> d = r1 - r2;
  AndoLab::Vector3d <double> dn = d.n();
  const double distance = d.abs();

  const double R_H = Radius_of_Earth + Altitude;

  const double hantei = (r2%dn)*(r2%dn) - r2.abs()*r2.abs() + R_H*R_H;

  /* 解がない */
  if ( hantei <= 0.0 ){
    num = 0;
  } else {

    double k_far = -1.0*r2%dn + std::sqrt(hantei);
    double k_near = -1.0*r2%dn - std::sqrt(hantei );

    /* r2->r1の領域に一つもAltitudeが入っていない */
    if ( (k_near >= distance) || (k_far <= 0.0)
        || ( (k_near <= 0.0) && (k_far >= distance) ) ){
      num = 0;

      /* r2->r1の領域に, Altitudeが1つ入っている */

      /* Nearが入っている */
    } else if ( (k_near > 0.0) && (k_near < distance) && (k_far > distance) ){
      num = 1;
      r[0] = r2 + k_near*dn;
    } else if ( (k_near < 0.0) && (k_far > 0.0) && (k_far < distance) ){
      num = 1;
      r[0] = r2 + k_far*dn;
    } else if ( (k_near > 0.0) && (k_far < distance) ){
      num = 2;
      r[0] = r2 + k_near*dn;
      r[1] = r2 + k_far*dn;
    } else {
      std::cout << "Error(Across): " << k_near << " & "
          << k_far << " ,, " << distance << std::endl;
      exit(0);
    }

  }
}

void Across_pts::copy(AndoLab::Vector3d <double> *dest){

  for(int i = 0; i < num; i++){
    dest[i] = r[i];
  }
}
